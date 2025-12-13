import random
import re
import json
import traceback
import base64
import asyncio
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple, Callable
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from langchain_openai import ChatOpenAI
from langchain_core.tools import Tool
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage

from src.tools.registry import TOOLS, auto_import_tools
from src.models.prompt import PromptTemplates

auto_import_tools()


class ConversationSession:
    """Manage a multi-turn conversation session."""
    def __init__(self, session_id: str, max_turns: int = 10):
        self.session_id = session_id
        self.max_turns = max_turns
        self.turn_count = 0
        self.history = []
        self.chat_history = []
        
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
        for i, turn in enumerate(self.history[-3:], 1):
            summary += f"Q{i}: {turn['user_query']}\n"
            if 'solution' in turn:
                solution_match = re.search(r"<solution>(.*?)</solution>", turn['solution'], re.S)
                solution_text = solution_match.group(1).strip() if solution_match else turn['solution']
                summary += f"A{i}: {solution_text}\n"
        return summary


class AgentWithTools:
    """
    Agent with real-time streaming support for frontend integration.
    """
    
    def __init__(self, tools: Dict[str, Dict], data_list: List[Dict[str, str]],
                 api_key: str, model_name: str = "deepseek-chat", temperature: float = 0.0, 
                 result_dir: str = "./results", max_turns_per_session: int = 10, 
                 max_code_iterations: int = 3):
        self.tools_meta = tools
        self.data_list = data_list
        self.llm = ChatOpenAI(model=model_name, temperature=temperature, api_key=api_key, 
                             base_url="https://api.deepseek.com")
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
        """Safely execute generated code in a sandboxed environment."""
        sandbox_globals = {"__builtins__": __builtins__}
        
        registered_count = 0
        for name, meta in self.tools_meta.items():
            if not isinstance(meta, dict):
                continue
                
            func = meta.get("func")
            if func is None or not callable(func):
                continue
            if not callable(func):
                continue

            sandbox_globals[name] = func
            registered_count += 1
        
        try:
            import matplotlib
            matplotlib.use('Agg')  # Use non-interactive backend
            import matplotlib.pyplot as plt
            import numpy as np
            sandbox_globals['plt'] = plt
            sandbox_globals['np'] = np
            plt.rcParams['font.sans-serif'] = ['Noto Sans CJK']
            plt.rcParams['axes.unicode_minus'] = False
        except ImportError:
            if code_type == "visualization":
                return {
                    "success": False,
                    "error": "matplotlib or numpy not available",
                    "error_type": "ImportError"
                }

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
                error_info["hint"] = "Function or variable not found."
            elif isinstance(e, TypeError):
                error_info["hint"] = "Type error. Check function arguments."
            elif isinstance(e, KeyError):
                error_info["hint"] = "Key not found in data structure."
            elif isinstance(e, IndexError):
                error_info["hint"] = "Index out of range. Check array dimensions."
            
            return error_info
    
    def _save_visualization(self, turn_id: str) -> Optional[str]:
        """Save generated matplotlib visualization."""
        try:
            if plt.get_fignums():
                vis_path = self.result_dir / f"visualization_{turn_id}.png"
                plt.savefig(vis_path, dpi=150, bbox_inches='tight')
                plt.close('all')
                return str(vis_path)
        except Exception as e:
            print(f"[Warning] Failed to save visualization: {e}")
        return None
    
    def _image_to_base64(self, image_path: str) -> Optional[str]:
        """Convert image to base64."""
        try:
            with open(image_path, 'rb') as f:
                return base64.b64encode(f.read()).decode('utf-8')
        except Exception as e:
            return None
    
    def _execute_with_retry(self, initial_code: str, context: dict, code_type: str,
                           full_specs: Dict[str, Any], user_query: str, language: str = 'en', status_callback: Optional[Callable] = None) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
        """
        Execute code with retry mechanism and status callbacks.
        """
        iteration_history = []
        current_code = initial_code

        for iteration in range(self.max_code_iterations):
            if status_callback:
                status_callback({
                    "type": "execution_attempt",
                    "code_type": code_type,
                    "iteration": iteration + 1,
                    "max_iterations": self.max_code_iterations
                })
            
            exec_result = self._safe_execute_tool(current_code, context, code_type)
            
            iteration_record = {
                "iteration": iteration + 1,
                "code": current_code,
                "success": exec_result["success"]
            }
            
            if exec_result["success"]:
                iteration_record["output"] = {k: str(v)[:] for k, v in exec_result["output"].items()}
                iteration_history.append(iteration_record)
                return exec_result, iteration_history
            
            iteration_record["error"] = exec_result["error"]
            iteration_record["error_type"] = exec_result["error_type"]
            iteration_history.append(iteration_record)
            
            if iteration < self.max_code_iterations - 1:
                if status_callback:
                    status_callback({
                        "type": "code_fix_attempt",
                        "code_type": code_type,
                        "error": exec_result["error"]
                    })
                
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
                else:
                    break
        
        return exec_result, iteration_history
    
    def create_session(self, session_id: Optional[str] = None) -> str:
        """Create a new conversation session."""
        if session_id is None:
            session_id = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        self.sessions[session_id] = ConversationSession(session_id, self.max_turns_per_session)
        return session_id
    
    def end_session(self, session_id: str):
        """End and cleanup session."""
        if session_id in self.sessions:
            del self.sessions[session_id]
    
    async def run_streaming(self, user_query: str, session_id: Optional[str] = None,
                           status_callback: Optional[Callable] = None) -> Dict[str, Any]:
        """
        Run agent with real-time status updates via callback.
        
        Args:
            - user_query: User's question
            - session_id: Optional session ID for conversation continuity
            - status_callback: Async function to send status updates to frontend
                           Format: async def callback(data: dict)
        """
        language = self._detect_language(user_query)

        # Create or get session
        if session_id is None:
            session_id = self.create_session()
        elif session_id not in self.sessions:
            session_id = self.create_session(session_id)
        
        session = self.sessions[session_id]
        
        if session.is_expired():
            self.end_session(session_id)
            session_id = self.create_session()
            session = self.sessions[session_id]
        
        result_record = {
            "session_id": session_id,
            "user_query": user_query,
            "language": language
        }
        
        # Status update helper
        async def send_status(data: dict):
            if status_callback:
                await status_callback(data)
        
        # Step 1: Planning
        await send_status({
            "stage": "planning",
            "message": "正在分析问题并选择数据和工具..."
        })
        
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

        system_prompt_1 = PromptTemplates.get_planning_prompt(
            language, tool_overview, data_overview, context_summary
        )
        planning_input = PromptTemplates.get_planning_input(language, user_query)

        planning_response = self.llm.invoke([
            SystemMessage(content=system_prompt_1),
            HumanMessage(content=planning_input)
        ])
        
        result_record["plan"] = planning_response.content

        # Extract selections
        selected_data_match = re.search(r"<selected_data>(.*?)</selected_data>", planning_response.content, re.S)
        selected_data_name = [d.strip() for d in selected_data_match.group(1).split(",")] if selected_data_match else []
        
        selected_data_path = []
        if selected_data_name:
            for name in selected_data_name:
                for data in self.data_list:
                    if data["name"] == name:
                        selected_data_path.append(data["path"])
        
        if not selected_data_path:
            selected_data_path = [self.data_list[0]["path"]]
            selected_data_name = [self.data_list[0]["name"]]
        
        selected_tool_match = re.search(r"<selected_tools>(.*?)</selected_tools>", planning_response.content, re.S)
        selected_tool_names = [t.strip() for t in selected_tool_match.group(1).split(",")] if selected_tool_match else []
        
        result_record["selected_data"] = {"name": selected_data_name, "path": selected_data_path}
        result_record["selected_tools"] = selected_tool_names

        # Send data and tool selection to frontend
        await send_status({
            "stage": "selection",
            "selected_data": selected_data_name,
            "selected_tools": selected_tool_names
        })

        # Step 2: Generate code
        await send_status({
            "stage": "code_generation",
            "message": "正在生成执行代码..."
        })
        
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
        
        result_record["generated_code"] = exec_response.content

        # Step 3: Execute computation
        await send_status({
            "stage": "execution",
            "message": "正在执行计算代码..."
        })
        
        compute_match = re.search(r"<compute>(.*?)</compute>", exec_response.content, re.S)
        if not compute_match:
            result_record["error"] = "No computation code found"
            session.add_turn(result_record)
            await send_status({
                "stage": "error",
                "error": "代码生成失败"
            })
            return result_record
        
        compute_code = compute_match.group(1).strip()
        turn_id = f"{session_id}_turn{session.turn_count + 1}"
        
        # Create callback wrapper for synchronous execution
        def sync_status_callback(data):
            asyncio.create_task(send_status(data))
        
        compute_result, compute_iterations = self._execute_with_retry(
            compute_code,
            context={"selected_data": selected_data_path[0] if selected_data_path else ""},
            code_type="compute",
            full_specs=full_specs,
            user_query=user_query,
            language=language,
            status_callback=sync_status_callback
        )

        result_record["compute_iterations"] = compute_iterations
        observation_parts = []
        
        if compute_result["success"]:
            compute_obs = {k: str(v)[:500] for k, v in compute_result["output"].items()}
            observation_parts.append({
                "type": "computation",
                "status": "success",
                "results": compute_obs
            })

            # Step 4: Visualization
            await send_status({
                "stage": "visualization",
                "message": "正在生成可视化..."
            })
            
            visualize_match = re.search(r"<visualize>(.*?)</visualize>", exec_response.content, re.S)
            if visualize_match:
                visualize_code = visualize_match.group(1).strip()
                
                if visualize_code:
                    vis_context = {"selected_data": selected_data_path[0] if selected_data_path else ""}
                    vis_context.update(compute_result["output"])
                    
                    vis_result, vis_iterations = self._execute_with_retry(
                        visualize_code,
                        context=vis_context,
                        code_type="visualization",
                        full_specs=full_specs,
                        user_query=user_query,
                        language=language,
                        status_callback=sync_status_callback
                    )
                    
                    result_record["visualization_iterations"] = vis_iterations
                    
                    if vis_result["success"]:
                        vis_path = self._save_visualization(turn_id)
                        if vis_path:
                            result_record["visualization_path"] = vis_path
                            vis_base64 = self._image_to_base64(vis_path)
                            result_record["visualization_base64"] = vis_base64
                            
                            # Send visualization to frontend
                            await send_status({
                                "stage": "visualization_ready",
                                "image_base64": vis_base64,
                                "image_path": vis_path
                            })
                            
                            observation_parts.append({
                                "type": "visualization",
                                "status": "success",
                                "path": vis_path
                            })
                        else:
                            observation_parts.append({
                                "type": "visualization",
                                "status": "not_applicable",
                                "message": "No visualization was generated for this query"
                            })
                    else:
                        observation_parts.append({
                            "type": "visualization",
                            "status": "failed",
                            "error": vis_result['error']
                        })
        else:
            observation_parts.append({
                "type": "computation",
                "status": "failed",
                "error": compute_result['error']
            })

        observation_text = json.dumps(observation_parts, indent=2, ensure_ascii=False)
        result_record["observation"] = observation_text

        # Step 5: Generate final answer
        await send_status({
            "stage": "generating_answer",
            "message": "正在生成最终回复..."
        })

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
        
        result_record["solution"] = final_response.content
        
        # Extract solution text
        solution_match = re.search(r"<solution>(.*?)</solution>", final_response.content, re.S)
        solution_text = solution_match.group(1).strip() if solution_match else final_response.content
        
        # Send final answer
        await send_status({
            "stage": "complete",
            "solution": solution_text
        })

        session.add_turn(result_record)
        return result_record