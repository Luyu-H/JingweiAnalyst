import random
import re
import json
import traceback
import base64
from datetime import datetime
from typing import Any, Dict, List, Optional
from pathlib import Path
import matplotlib.pyplot as plt

from langchain_openai import ChatOpenAI
from langchain_core.tools import Tool
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage
from langchain_core.prompts import PromptTemplate

from src.tools.registry import TOOLS, auto_import_tools, export_tools_doc

# import all tools and register them to TOOLS
auto_import_tools()


class ConversationSession:
    """Manage a multi-turn conversation session."""
    def __init__(self, session_id: str, max_turns: int = 10):
        self.session_id = session_id
        self.max_turns = max_turns
        self.turn_count = 0
        self.history = []  # Record of all turns
        self.chat_history = []  # Record for LLM context
        
    def add_turn(self, turn_record: Dict[str, Any]):
        """Add a turn record to the session history."""
        self.turn_count += 1
        turn_record["turn_number"] = self.turn_count
        turn_record["timestamp"] = datetime.now().isoformat()
        self.history.append(turn_record)
        
        self.chat_history.append(HumanMessage(content=turn_record["user_query"]))
        if "solution" in turn_record:
            self.chat_history.append(AIMessage(content=turn_record["solution"]))
    
    def is_expired(self) -> bool:
        return self.turn_count >= self.max_turns
    
    def get_context_summary(self) -> str:
        """Generate a brief summary of the conversation history for context."""
        if not self.history:
            return ""
        
        summary = f"Previous conversation (Turn {self.turn_count}):\n"
        for i, turn in enumerate(self.history[-3:], 1):  # only summarize last 3 turns
            summary += f"Q{i}: {turn['user_query']}\n"
            if 'solution' in turn:
                solution_match = re.search(r"<solution>(.*?)</solution>", turn['solution'], re.S)
                solution_text = solution_match.group(1).strip() if solution_match else turn['solution']
                # summary += f"A{i}: {solution_text[:200]}...\n"  # limit length
                summary += f"A{i}: {solution_text}\n"
        return summary


class AgentWithTools:
    """An agent that can select and call provided tools to answer scientific questions with multi-turn conversation support."""
    
    def __init__(self, tools: Dict[str, Dict], data_list: List[Dict[str, str]],
                 api_key: str, model_name: str = "deepseek-chat", temperature: float = 0.0, 
                 result_dir: str = "./results", max_turns_per_session: int = 10, max_code_iterations: int = 3):
        self.tools_meta = tools
        self.data_list = data_list
        self.llm = ChatOpenAI(model=model_name, temperature=temperature, api_key=api_key, base_url="https://api.deepseek.com")
        self.tools = self._build_tool_list(tools)
        self.result_dir = Path(result_dir)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        
        self.sessions: Dict[str, ConversationSession] = {}
        self.max_turns_per_session = max_turns_per_session
    
    def _build_tool_list(self, tools_meta: Dict[str, Dict]) -> List[Tool]:
        tool_list = []
        for name, meta in tools_meta.items():
            tool = Tool.from_function(
                func=meta["func"],
                name=name,
                description=meta.get("description", ""),
                return_direct=False
            )
            tool_list.append(tool)
        return tool_list

    def _safe_execute_tool(self, code: str, context: dict = None) -> Dict[str, Any]:
        """Safely execute generated code in a sandboxed environment with registered tools."""
        sandbox_globals = {"__builtins__": __builtins__}
        
        registered_count = 0
        for name, meta in self.tools_meta.items():
            if not isinstance(meta, dict):
                print(f"[Warning] Tool '{name}' meta is not a dict: {type(meta)}")
                continue
                
            func = meta.get("func")
            if func is None:
                print(f"[Warning] Tool '{name}' has no 'func' key in meta")
                continue
            if not callable(func):
                print(f"[Warning] Tool '{name}' func is not callable: {type(func)}")
                continue
            
            sandbox_globals[name] = func
            registered_count += 1
        
        # Add matplotlib and numpy for visualization
        try:
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import matplotlib.pyplot as plt
            import numpy as np
            sandbox_globals['plt'] = plt
            sandbox_globals['np'] = np
        except ImportError:
            print("[Warning] matplotlib or numpy not available for visualization")
        
        print(f"[Info] Successfully registered {registered_count} tools to sandbox")

        if context:
            sandbox_globals.update(context)
            
        sandbox_locals = {}
        try:
            exec(code, sandbox_globals, sandbox_locals)
            return {"success": True, "output": sandbox_locals}
        except Exception as e:
            tb = traceback.format_exc()
            error_info = {
                "success": False,
                "error": str(e),
                "error_type": type(e).__name__,
                "traceback": tb,
                "registered_tools": list(self.tools_meta.keys())
            }
            
            if isinstance(e, NameError):
                error_info["hint"] = "Function not found. Check if the function name is correct and registered."
            if isinstance(e, TypeError) and "'dict' object is not callable" in str(e):
                error_info["hint"] = "Trying to call a dict as function. Check function registration."
            
            return error_info
    
    def _save_visualization(self, locals_dict: Dict[str, Any], turn_id: str) -> Optional[str]:
        """Save generated matplotlib visualization if exists."""
        try:
            if plt.get_fignums():
                vis_path = self.result_dir / f"visualization_{turn_id}.png"
                plt.savefig(vis_path, dpi=150, bbox_inches='tight')
                plt.close('all')
                print(f"[System] Visualization saved to {vis_path}")
                return str(vis_path)
        except Exception as e:
            print(f"[Warning] Failed to save visualization: {e}")
        
        return None
    
    def _image_to_base64(self, image_path: str) -> Optional[str]:
        """Convert an image file to base64 string."""
        try:
            with open(image_path, 'rb') as f:
                return base64.b64encode(f.read()).decode('utf-8')
        except Exception as e:
            print(f"[Warning] Failed to encode image: {e}")
            return None
    
    def _save_session(self, session: ConversationSession):
        """Save session history to a JSON file."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = self.result_dir / f"session_{session.session_id}_{timestamp}.json"
        
        session_data = {
            "session_id": session.session_id,
            "total_turns": session.turn_count,
            "max_turns": session.max_turns,
            "history": session.history
        }
        
        with open(fname, "w", encoding="utf-8") as f:
            json.dump(session_data, f, indent=2, ensure_ascii=False)
        print(f"[System] Session saved to {fname}")
    
    def create_session(self, session_id: Optional[str] = None) -> str:
        """Create a new conversation session."""
        if session_id is None:
            session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        self.sessions[session_id] = ConversationSession(session_id, self.max_turns_per_session)
        print(f"[System] Created new session: {session_id}")
        return session_id
    
    def end_session(self, session_id: str):
        """End a conversation session and save its history."""
        if session_id in self.sessions:
            self._save_session(self.sessions[session_id])
            del self.sessions[session_id]
            print(f"[System] Session {session_id} ended and saved")
    
    def run(self, user_query: str, session_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Run the agent to answer a user query within a conversation session.
        
        Args:
            - user_query, str: User's question
            - session_id: Optional session ID for multi-turn conversation
            
        Returns:
            Dict, the result record including plan, selected data/tools, generated code, observation, and solution.
        """
        if session_id is None:
            session_id = self.create_session()
        elif session_id not in self.sessions:
            print(f"[Warning] Session {session_id} not found, creating new session")
            session_id = self.create_session(session_id)
        
        session = self.sessions[session_id]
        
        if session.is_expired():
            print(f"[Warning] Session {session_id} has reached max turns, creating new session")
            self.end_session(session_id)
            session_id = self.create_session()
            session = self.sessions[session_id]
        
        result_record = {
            "session_id": session_id,
            "user_query": user_query
        }
        
        context_summary = session.get_context_summary()
        
        tool_overview = {
            k: {"description": v["description"]}
            for k, v in self.tools_meta.items()
        }
        
        data_overview = {
            data["name"]: {
                "description": data.get("description", ""),
                "path": data["path"]
            }
            for data in self.data_list
        }

        # Step 1: Planning, including selecting data and tools
        system_prompt_1 = f"""You are a scientific data analysis assistant with conversation history support.

{context_summary}

Available tools:
{json.dumps(tool_overview, indent=2)}

Available datasets:
{json.dumps(data_overview, indent=2)}
Note that if the user does not specify which version of the data is needed, select the latest version of relevant variables.

Analyze the user's question and select appropriate data and tools. Follow this format strictly:
<plan>step-by-step reasoning about which data and tools are needed</plan>
<selected_data>comma-separated list of dataset names you will use</selected_data>
<selected_tools>comma-separated list of tool names you will use</selected_tools>
"""

        planning_input = f"User question: {user_query}\n\nNow output your reasoning and selections."
        planning_response = self.llm.invoke([
            SystemMessage(content=system_prompt_1),
            HumanMessage(content=planning_input)
        ])
        print("\n=== Plan & Selection ===\n", planning_response.content)
        result_record["plan"] = planning_response.content

        # Extract selected data
        selected_data_match = re.search(r"<selected_data>(.*?)</selected_data>", planning_response.content, re.S)
        selected_data_name = [d.strip() for d in selected_data_match.group(1).split(",")] if selected_data_match else []
        
        selected_data_path = []
        if selected_data_name:
            for name in selected_data_name:
                for data in self.data_list:
                    if data["name"] == name:
                        selected_data_path.append(data["path"])
        
        if not selected_data_path:
            print(f"[Warning] Could not find data '{selected_data_name}', using first dataset")
            selected_data_path = self.data_list[0]["path"]
            selected_data_name = self.data_list[0]["name"]
        
        print(f"[System] Selected data: {selected_data_name} ({selected_data_path})")
        result_record["selected_data"] = {"name": selected_data_name, "path": selected_data_path}

        selected_tool_match = re.search(r"<selected_tools>(.*?)</selected_tools>", planning_response.content, re.S)
        selected_tool_names = [t.strip() for t in selected_tool_match.group(1).split(",")] if selected_tool_match else []
        result_record["selected_tools"] = selected_tool_names

        # Step 2: Generate executable code
        full_specs = {
            name: {
                "description": self.tools_meta[name]["description"],
                "args": self.tools_meta[name].get("args", {}),
                "returns": self.tools_meta[name].get("returns", "")
            }
            for name in selected_tool_names if name in self.tools_meta
        }

        system_prompt_2 = f"""Here are the full specifications for the tools you have selected:

{json.dumps(full_specs, indent=2)}

Generate executable Python code between <execute></execute> tags to solve: "{user_query}"

Requirements:
1. The dataset path is available as variable `selected_data` = {str(selected_data_path)}"
2. Call the necessary tools to perform calculations
3. Store key results in variables (e.g., result, mean_value, etc.)
4. If the result can be visualized (time series, spatial distribution, histogram, etc.), 
   generate matplotlib code to create an appropriate visualization
5. Available libraries: matplotlib.pyplot as plt, numpy as np
6. Use plt.savefig() is NOT needed - the system will handle saving automatically
7. Add proper and scientific labels, title, and legend to your plots

Example structure:
<execute>
# Step 1: Load and process data
result = some_tool(selected_data, ...)

# Step 2: Perform calculations
mean_value = result.mean()

# Step 3: Create visualization (if appropriate)
plt.figure(figsize=(10, 6))
plt.plot(result)
plt.xlabel('...')
plt.ylabel('...')
plt.title('...')
plt.grid(True)
</execute>
"""
        
        exec_response = self.llm.invoke([
            SystemMessage(content=system_prompt_2),
            HumanMessage(content="Generate code now.")
        ])
        print("\n=== Generated Code ===\n", exec_response.content)
        result_record["generated_code"] = exec_response.content

        # Step 3: Extract and execute code
        code_match = re.search(r"<execute>(.*?)</execute>", exec_response.content, re.S)
        if not code_match:
            print("[Error] No <execute> code block found.")
            result_record["error"] = "No executable code found"
            session.add_turn(result_record)
            return result_record
        
        code_block = code_match.group(1).strip()
        turn_id = f"{session_id}_turn{session.turn_count + 1}"
        exec_result = self._safe_execute_tool(code_block, context={"selected_data": selected_data_path})

        if exec_result["success"]:
            print("\n=== Observation ===")
            obs = {k: str(v)[:] for k, v in exec_result["output"].items()}
            print(json.dumps(obs, indent=2))
            observation_text = json.dumps(obs, indent=2)
            
            vis_path = self._save_visualization(exec_result["output"], turn_id)
            if vis_path:
                result_record["visualization_path"] = vis_path
                result_record["visualization_base64"] = self._image_to_base64(vis_path)
        else:
            print("\n=== Execution Error ===")
            print(exec_result["traceback"])
            observation_text = f"Execution failed: {exec_result['error']}\n{exec_result.get('hint', '')}"

        result_record["observation"] = observation_text

        # Step 4: Generate final answer
        solution_prompt = f"""
Here is the result of your code execution:

<observation>
{observation_text}
</observation>

{"A visualization has been generated and saved." if "visualization_path" in result_record else ""}

Please provide a clear scientific summary inside <solution></solution> tags.
Include:
1. Direct answer to the question
2. Key numerical results
3. Description of visualization (if generated)
4. Any important insights or patterns observed
"""
        final_response = self.llm.invoke([
            SystemMessage(content=solution_prompt),
            HumanMessage(content="Provide the final explanation.")
        ])
        print("\n=== Final Solution ===\n", final_response.content)
        result_record["solution"] = final_response.content

        # Update session history
        session.add_turn(result_record)
        
        return result_record
    

if __name__ == "__main__":
    data_list = [
        {
            "name": "nitrate_data_v1.1",
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_no3_sjtu_jingwei_v1_1.zarr",
            "description": "Nitrate concentration data from ocean reconstruction, version 1.1"
        },
        {
            "name": "oxygen_data_v1.1",
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr",
            "description": "Dissolved oxygen concentration data from ocean reconstruction, version 1.1"
        },
        {
            "name": "phosphate_data_v1.1", 
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr",
            "description": "Phosphate concentration data from ocean reconstruction, version 1.1"
        },
        {
            "name": "silicate_data_v1.1", 
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_si_sjtu_jingwei_v1_1.zarr",
            "description": "Silicate concentration data from ocean reconstruction, version 1.1"
        },
        {
            "name": "temperature_data_v1.1", 
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_temp_sjtu_jingwei_v1_1.zarr",
            "description": "Sea temperature data from ocean reconstruction, version 1.1"
        },
        {
            "name": "salinity_data_v1.1", 
            "path": "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_so_sjtu_jingwei_v1_1.zarr",
            "description": "Sea salinity data from ocean reconstruction, version 1.1"
        }
    ]

    api_key = "sk-edc62146735146e4a9904740057ace47"
    agent = AgentWithTools(TOOLS, data_list, api_key, max_turns_per_session=10)
    
    session_id = agent.create_session()
    
    # first question
    result1 = agent.run(
        "What is the mean oxygen concentration between 2000 and 2010 at surface level?",
        session_id=session_id
    )
    # second question in the same session
    result2 = agent.run(
        "How does this compare to the phosphate levels in the same period?",
        session_id=session_id
    )
    
    agent.end_session(session_id)