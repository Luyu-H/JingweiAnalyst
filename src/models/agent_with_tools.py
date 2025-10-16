import random
import re
import json
import traceback
from datetime import datetime
from typing import Any, Dict, List

from langchain_openai import ChatOpenAI
from langchain.agents import Tool, initialize_agent
from langchain.schema import SystemMessage, HumanMessage
from langchain.prompts import PromptTemplate
from langchain.memory import ConversationBufferMemory

from src.tools.registry import TOOLS, auto_import_tools, export_tools_doc

# import all tools and register them to TOOLS
auto_import_tools()


class AgentWithTools:
    """An agent that can select and call provided tools to answer scitific questions."""
    def __init__(self, tools: Dict[str, Dict], api_key: str, model_name: str = "deepseek-chat", temperature: float = 0.0, result_dir: str = "./results"):
        self.tools_meta = tools
        self.llm = ChatOpenAI(
            model=model_name,
            temperature=temperature,
            api_key=api_key,
            base_url="https://api.deepseek.com"
        )
        self.tools = self._build_tool_list(tools)
        self.memory = ConversationBufferMemory(memory_key="chat_history")
        self.result_dir = result_dir
    
    def _select_data(self, data_list:List[str]) -> str:
        # return random.choice(data_list)
        return data_list[2]
    
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
            
            # 成功注册函数
            sandbox_globals[name] = func
            registered_count += 1
        
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
                error_info["hint"] = "Trying to call a dict as function. Possible causes: 1) Function not properly registered, 2) Variable name conflicts with dict."
            
            return error_info
    
    def _save_result(self, result_dict: Dict[str, Any], prefix="agent_output"):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = f"{self.result_dir}/{prefix}_{timestamp}.json"
        with open(fname, "w", encoding="utf-8") as f:
            json.dump(result_dict, f, indent=2, ensure_ascii=False)
        print(f"[System] Result saved to {fname}")
    
    def run(self, user_query: str, data_list: List[str]) -> str:
        result_record = {"user_query": user_query}

        selected_data = self._select_data(data_list)
        print(f"[System] Selected data table: {selected_data}")
        result_record["selected_data"] = selected_data

        tool_overview = {
            k: {"description": v["description"]}
            for k, v in self.tools_meta.items()
        }

        # Step 1: Construct system prompt for planning and tool selection
        system_prompt_1 = f"""You are a scientific data analysis assistant.
You have access to the following tools (only names and descriptions shown below):

{json.dumps(tool_overview, indent=2)}

Follow this format strictly:
<plan> step-by-step reasoning </plan>
<selected_tools> list the names of tools you will use </selected_tools>
"""

        plan_prompt = PromptTemplate(
            input_variables=["question"],
            template="{system_prompt}\nUser question:\n{question}\nNow output your reasoning."
        )

        planning_input = plan_prompt.format(system_prompt=system_prompt_1, question=user_query)
        planning_response = self.llm.invoke([
            SystemMessage(content=system_prompt_1),
            HumanMessage(content=planning_input)
        ])
        print("\n=== Plan & Tool Selection ===\n", planning_response.content)
        result_record["plan"] = planning_response.content

        selected_tool_names = re.findall(r"<selected_tools>(.*?)</selected_tools>", planning_response.content, re.S)
        selected_tool_names = [t.strip() for t in selected_tool_names[0].split(",")] if selected_tool_names else []

        # Step 2: Construct system prompt for tool execution
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

Now generate executable Python code between <execute></execute> tags to solve the question "{user_query}".
You can directly call these tools by their names (functions are already available in the environment).
The dataset path has been provided as variable `selected_data`, and you can use it directly without definition.
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
            self._save_result(result_record)
            return result_record
        
        code_block = code_match.group(1).strip()
        # code_block = f'selected_data = "{selected_data}"\n' + code_block
        exec_result = self._safe_execute_tool(code_block, context={"selected_data": selected_data})

        if exec_result["success"]:
            print("\n=== Observation ===")
            obs = {k: str(v) for k, v in exec_result["output"].items()}
            print(json.dumps(obs, indent=2))
            observation_text = json.dumps(obs, indent=2)
        else:
            print("\n=== Execution Error ===")
            print(exec_result["traceback"])
            observation_text = f"Execution failed: {exec_result['error']}"

        result_record["observation"] = observation_text

        # Step 4: Generate final answer
        solution_prompt = f"""
Here is the result of your code execution:

<observation>
{observation_text}
</observation>

Please summarize the scientific result clearly inside <solution></solution>.
"""
        final_response = self.llm.invoke([
            SystemMessage(content=solution_prompt),
            HumanMessage(content="Provide the final explanation.")
        ])
        print("\n=== Final Solution ===\n", final_response.content)
        result_record["solution"] = final_response.content

        # Save all results
        self._save_result(result_record)
        return result_record
    

if __name__ == "__main__":
    data_list = [
        "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_no3_sjtu_jingwei_v1_1.zarr",
        "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_po4_sjtu_jingwei_v1_1.zarr",
        "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"
    ]

    api_key = "sk-edc62146735146e4a9904740057ace47"
    agent = AgentWithTools(TOOLS, api_key)
    query = "What is the mean oxygen concentration between 2000 and 2010 at surface level?"
    result = agent.run(query, data_list)