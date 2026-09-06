import random
import re
import os
import json
import traceback
import base64
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple
from pathlib import Path
import matplotlib.pyplot as plt

from langchain_openai import ChatOpenAI
from langchain_core.tools import Tool
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage

from src.tools.registry import TOOLS, auto_import_tools, export_tools_doc
from src.models.prompt import PromptTemplates

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
                 api_key: str, model_name: str = "deepseek-reasoner", temperature: float = 0.0, 
                 result_dir: str = "./results", max_turns_per_session: int = 10, max_code_iterations: int = 3):
        self.tools_meta = tools
        self.data_list = data_list
        self.llm = ChatOpenAI(model=model_name, temperature=temperature, api_key=api_key, base_url="https://api.deepseek.com")
        # self.llm = ChatOpenAI(model=model_name, temperature=temperature, api_key=api_key)
        self.tools = self._build_tool_list(tools)
        self.result_dir = Path(result_dir)
        self.result_dir.mkdir(parents=True, exist_ok=True)
        
        self.sessions: Dict[str, ConversationSession] = {}
        self.max_turns_per_session = max_turns_per_session
        self.max_code_iterations = max_code_iterations
    
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

    def _detect_language(self, text: str) -> str:
        """Detect if the input text is Chinese or English.
        Returns 'zh' for Chinese, 'en' for English."""
        chinese_char_count = len(re.findall(r'[\u4e00-\u9fff]', text))
        english_char_count = len(re.findall(r'[a-zA-Z]', text))

        if chinese_char_count > english_char_count:
            return 'zh'
        else:
            return 'en'

    def _safe_execute_tool(self, code: str, context: dict = None, code_type: str = "compute") -> Dict[str, Any]:
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
            plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC']
            plt.rcParams['axes.unicode_minus'] = False
            
        except ImportError:
            if code_type == "visualization":
                return {
                    "success": False,
                    "error": "matplotlib or numpy not available",
                    "error_type": "ImportError"
                }
        
        print(f"[Info] Successfully registered {registered_count} tools to sandbox")

        if code_type == "compute":
            print(f"[Info] Executing computation code with {registered_count} registered tools")
        else:
            print(f"[Info] Executing visualization code")

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
                error_info["hint"] = "Function or variable not found. Check spelling and availability."
            elif isinstance(e, TypeError):
                error_info["hint"] = "Type error. Check function arguments and data types."
            elif isinstance(e, KeyError):
                error_info["hint"] = "Key not found in data structure. Check data keys."
            elif isinstance(e, IndexError):
                error_info["hint"] = "Index out of range. Check array dimensions."
            
            return error_info
    
    def _save_visualization(self, turn_id: str) -> Optional[str]:
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
    
    def _execute_with_retry(self, initial_code: str, context: dict, code_type: str,
                           full_specs: Dict[str, Any], user_query: str, language: str = 'en') -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        """
        Execute code with retries upon failure, regenerating code as needed.
        
        Returns:
            - Dict: Final execution result dict
            - List: List of all execution attempts with code and results
        """
        iteration_history = []
        current_code = initial_code

        for iteration in range(self.max_code_iterations):
            print(f"\n--- Iteration {iteration + 1}/{self.max_code_iterations} ({code_type}) ---")
            
            exec_result = self._safe_execute_tool(current_code, context, code_type)
            
            iteration_record = {
                "iteration": iteration + 1,
                "code": current_code,
                "success": exec_result["success"]
            }
            
            if exec_result["success"]:
                iteration_record["output"] = {k: str(v)[:] for k, v in exec_result["output"].items()}
                iteration_history.append(iteration_record)
                print(f"[Success] {code_type.capitalize()} code executed successfully!")
                return exec_result, iteration_history
            
            # If failed, record error info
            iteration_record["error"] = exec_result["error"]
            iteration_record["error_type"] = exec_result["error_type"]
            iteration_history.append(iteration_record)
            
            print(f"[Error] Execution failed: {exec_result['error']}")
            
            # Request code fix from LLM
            if iteration < self.max_code_iterations - 1:
                print("[System] Requesting code fix from LLM...")
                
                fix_prompt = PromptTemplates.get_fix_prompt(
                    language, current_code, exec_result, user_query, full_specs, code_type
                )
                system_content = PromptTemplates.get_fix_system_message(language)
                
                fix_response = self.llm.invoke([
                    SystemMessage(content=system_content),
                    HumanMessage(content=fix_prompt)
                ])
                
                code_match = re.search(r"<execute>(.*?)</execute>", fix_response.content, re.S)
                if code_match:
                    current_code = code_match.group(1).strip()
                    print(f"[System] Code fixed, retrying...")
                else:
                    print("[Warning] LLM did not return fixed code in proper format")
                    break
            else:
                print(f"[System] Max iterations reached for {code_type}")
        
        return exec_result, iteration_history
    
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
        language = self._detect_language(user_query)
        print(f"[System] Detected language: {'Chinese' if language == 'zh' else 'English'}")

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
            "user_query": user_query,
            "language": language
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
        system_prompt_1 = PromptTemplates.get_planning_prompt(
            language, tool_overview, data_overview, context_summary
        )
        planning_input = PromptTemplates.get_planning_input(language, user_query)

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

        system_prompt_2 = PromptTemplates.get_code_generation_prompt(
            language, full_specs, user_query, selected_data_path
        )
        exec_input = PromptTemplates.get_code_generation_input(language)
        
        exec_response = self.llm.invoke([
            SystemMessage(content=system_prompt_2),
            HumanMessage(content=exec_input)
        ])
        print("\n=== Generated Code ===\n", exec_response.content)
        result_record["generated_code"] = exec_response.content

        # Step 3: Extract and execute the computation code
        compute_match = re.search(r"<compute>(.*?)</compute>", exec_response.content, re.S)
        if not compute_match:
            print("[Error] No <compute> code block found.")
            result_record["error"] = "No computation code found"
            session.add_turn(result_record)
            return result_record
        
        compute_code = compute_match.group(1).strip()
        turn_id = f"{session_id}_turn{session.turn_count + 1}"
        
        compute_result, compute_iterations = self._execute_with_retry(
            compute_code,
            context={"selected_data": selected_data_path},
            code_type="compute",
            full_specs=full_specs,
            user_query=user_query,
            language=language
        )

        result_record["compute_iterations"] = compute_iterations

        observation_parts = []
        
        if compute_result["success"]:
            print("\n=== Computation Observation ===")
            compute_obs = {k: str(v)[:500] for k, v in compute_result["output"].items()}
            print(json.dumps(compute_obs, indent=2))
            observation_parts.append({
                "type": "computation",
                "status": "success",
                "results": compute_obs
            })

            # Step 4: Extract and execute the visualization code (if any)
            visualize_match = re.search(r"<visualize>(.*?)</visualize>", exec_response.content, re.S)
            if visualize_match:
                visualize_code = visualize_match.group(1).strip()
                
                if visualize_code:
                    vis_context = {"selected_data": selected_data_path}
                    vis_context.update(compute_result["output"])
                    
                    vis_result, vis_iterations = self._execute_with_retry(
                        visualize_code,
                        context=vis_context,
                        code_type="visualization",
                        full_specs=full_specs,
                        user_query=user_query,
                        language=language
                    )
                    
                    result_record["visualization_iterations"] = vis_iterations
                    
                    if vis_result["success"]:
                        print("\n=== Visualization Observation ===")
                        print("Visualization code executed successfully")
                        
                        vis_path = self._save_visualization(turn_id)
                        if vis_path:
                            result_record["visualization_path"] = vis_path
                            result_record["visualization_base64"] = self._image_to_base64(vis_path)
                            observation_parts.append({
                                "type": "visualization",
                                "status": "success",
                                "path": vis_path
                            })
                        else:
                            observation_parts.append({
                                "type": "visualization",
                                "status": "no_figure",
                                "message": "Code executed but no figure was generated"
                            })
                    else:
                        print("\n=== Visualization Failed ===")
                        print(f"Error: {vis_result['error']}")
                        observation_parts.append({
                            "type": "visualization",
                            "status": "failed",
                            "error": vis_result['error'],
                            "message": "Visualization failed but computation succeeded"
                        })
                else:
                    print("[Info] No visualization code provided")
                    observation_parts.append({
                        "type": "visualization",
                        "status": "not_applicable",
                        "message": "No visualization was generated for this query"
                    })
            else:
                observation_parts.append({
                    "type": "visualization",
                    "status": "not_provided",
                    "message": "No visualization block found in generated code"
                })
        else:
            print("\n=== Computation Failed ===")
            print(compute_result["traceback"])
            observation_parts.append({
                "type": "computation",
                "status": "failed",
                "error": compute_result['error'],
                "iterations": len(compute_iterations)
            })

        observation_text = json.dumps(observation_parts, indent=2, ensure_ascii=False)
        result_record["observation"] = observation_text

        # Step 4: Generate final answer
        has_computation_success = observation_parts[0]["status"] == "success"
        has_visualization = any(p.get("type") == "visualization" and p.get("status") == "success" 
                               for p in observation_parts)
        
        solution_prompt = PromptTemplates.get_solution_prompt(
            language, observation_text, has_computation_success, has_visualization
        )
        final_input = PromptTemplates.get_solution_input(language)
        
        final_response = self.llm.invoke([
            SystemMessage(content=solution_prompt),
            HumanMessage(content=final_input)
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

    def run_tests_from_folder(folder_path, agent):
        """
        Read all JSON files inside folder_path.
        Each JSON file contains an array of test items.
        For each item: create a new agent session and run the query.
        """

        # list all .json files
        json_files = [
            os.path.join(folder_path, f)
            for f in os.listdir(folder_path)
            if f.endswith(".json")
        ]

        print(f"Found {len(json_files)} test files.")

        all_results = []   # collect everything

        for json_file in json_files:
            print(f"\n=== Loading test file: {json_file} ===")
            with open(json_file, "r") as f:
                test_items = json.load(f)

            # test_items should be a list
            if not isinstance(test_items, list):
                raise ValueError(f"JSON file {json_file} must contain a list/array.")

            for idx, item in enumerate(test_items):
                query = item.get("query")

                if query is None:
                    print(f"Skipping item {idx} because no 'query' field.")
                    continue

                print(f"\n--- Running test case {idx+1}/{len(test_items)} ---")
                print(f"Query: {query}")

                # Create a new session for EACH query
                session_id = agent.create_session()

                try:
                    result = agent.run(query, session_id=session_id)
                    all_results.append({
                        "file": json_file,
                        "index": idx,
                        "query": query,
                        "result": result
                    })
                except Exception as e:
                    print(f"Error running test case {idx}: {e}")
                finally:
                    agent.end_session(session_id)

        return all_results

    api_key = "sk-edc62146735146e4a9904740057ace47"
    agent = AgentWithTools(TOOLS, data_list, api_key, model_name="deepseek-reasoner",
                           max_turns_per_session=10, max_code_iterations=3)
    
    folder = "/home/lorraine/jingwei_QA/data/web_applications"
    for _ in range(10):
        results = run_tests_from_folder(folder, agent)

    print("Finished all tests.")
    
    # session_id = agent.create_session()
    
    # # first question
    # result1 = agent.run(
    #     "In the North Pacific, evaluate the lagged correlation between mixed layer depth (MLD) and surface (10 m) dissolved oxygen anomalies on a seasonal scale for the period 1990–2010 (to identify lead–lag relationships), and construct a multivariate regression model to examine the contributions of MLD, sea-surface temperature (SST), and salinity to surface oxygen. Pay attention to checking the data dimensions. The longitude and latitude, depth, and time resolution of the initial data of all variables are the same.",
    #     session_id=session_id
    # )
    # # second question in the same session
    # result2 = agent.run(
    #     "How does this compare to the phosphate levels in the same period?",
    #     session_id=session_id
    # )
    
    # agent.end_session(session_id)