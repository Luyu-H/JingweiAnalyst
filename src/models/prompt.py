import json
from typing import Dict, Any


class PromptTemplates:
    """Collection of all prompt templates for the agent."""
    
    @staticmethod
    def get_planning_prompt(language: str, tool_overview: Dict, data_overview: Dict, context_summary: str) -> str:

        if language == 'zh':
            return f"""你是一个科学数据分析助手，支持对话历史功能。

{context_summary}

可用工具：
{json.dumps(tool_overview, indent=2, ensure_ascii=False)}

可用数据集：
{json.dumps(data_overview, indent=2, ensure_ascii=False)}
注意：如果用户没有指定需要哪个版本的数据，请选择相关变量的最新版本。

分析用户的问题并选择合适的数据和工具。严格按照以下格式输出：
<plan>关于需要哪些数据和工具的逐步推理</plan>
<selected_data>您将使用的数据集名称，用逗号分隔</selected_data>
<selected_tools>您将使用的工具名称，用逗号分隔</selected_tools>
"""
        
        else:  # English
            return f"""You are a scientific data analysis assistant with conversation history support.

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
    
    @staticmethod
    def get_planning_input(language: str, user_query: str) -> str:
        if language == 'zh':
            return f"用户问题：{user_query}\n\n现在输出您的推理和选择。"
        else:
            return f"User question: {user_query}\n\nNow output your reasoning and selections."
    
    @staticmethod
    def get_code_generation_prompt(language: str, full_specs: Dict, user_query: str, selected_data_path: str) -> str:

        if language == 'zh':
            return f"""以下是您选择的工具的完整规格说明：

{json.dumps(full_specs, indent=2, ensure_ascii=False)}

生成两个独立的代码块来解决："{user_query}"

要求：
1. 数据集路径可通过变量 `selected_data` = "{selected_data_path}" 获取
2. 可用库：matplotlib.pyplot as plt, numpy as np

按以下格式生成代码：

<compute>
# 仅计算代码
# 调用工具加载和处理数据
# 执行计算并将结果存储在变量中
# 示例：result, mean_value, trend 等
</compute>

<visualize>
# 仅可视化代码（此部分可选）
# 使用计算结果（它们将作为变量可用）
# 创建带有适当标签和标题的matplotlib图表
# 使用 plt.figure(), plt.plot(), plt.xlabel() 等
# 不要使用 plt.savefig() - 会自动处理
</visualize>

重要提示：
- 计算代码必须完整且独立
- 可视化代码可以引用计算中的变量
- 如果此查询不适合可视化，请将 <visualize> 留空或省略
"""
        
        else:  # English
            return f"""Here are the full specifications for the tools you have selected:

{json.dumps(full_specs, indent=2)}

Generate TWO separate code blocks to solve: "{user_query}"

Requirements:
1. The dataset path is available as variable `selected_data` = "{selected_data_path}"
2. Available libraries: matplotlib.pyplot as plt, numpy as np

Generate code in this format:

<compute>
# Computation code only
# Call tools to load and process data
# Perform calculations and store results in variables
# Examples: result, mean_value, trend, etc.
</compute>

<visualize>
# Visualization code only (this section is optional)
# Use results from computation (they will be available as variables)
# Create matplotlib plots with proper labels and titles
# Use plt.figure(), plt.plot(), plt.xlabel(), etc.
# DO NOT use plt.savefig() - handled automatically
</visualize>

Important:
- Computation code MUST be complete and standalone
- Visualization code can reference variables from computation
- If visualization is not appropriate for this query, leave <visualize> empty or omit it
"""
    
    @staticmethod
    def get_code_generation_input(language: str) -> str:
        if language == 'zh':
            return "现在生成两个代码块。"
        else:
            return "Generate the two code blocks now."
    
    @staticmethod
    def get_fix_prompt(language: str, code: str, exec_result: Dict, user_query: str, full_specs: Dict, code_type: str) -> str:

        if language == 'zh':
            return f"""以下{code_type}代码执行失败：

<failed_code>
{code}
</failed_code>

<error_info>
错误类型：{exec_result['error_type']}
错误信息：{exec_result['error']}
{exec_result.get('hint', '')}

追踪信息：
{exec_result.get('traceback', '不适用')}
</error_info>

原始任务：{user_query}

可用工具及其规格说明：
{json.dumps(full_specs, indent=2, ensure_ascii=False)}

请修复代码并在 <execute></execute> 标签内输出修正后的版本。
重点关注：
1. 修复上述具体错误
2. 确保正确的变量名和函数调用
3. 使用正确的语法和库函数
"""
        
        else:  # English
            return f"""The following {code_type} code failed to execute:

<failed_code>
{code}
</failed_code>

<error_info>
Error Type: {exec_result['error_type']}
Error Message: {exec_result['error']}
{exec_result.get('hint', '')}

Traceback:
{exec_result.get('traceback', 'N/A')}
</error_info>

Original task: {user_query}

Available tools and their specifications:
{json.dumps(full_specs, indent=2)}

Please fix the code and output the corrected version inside <execute></execute> tags.
Focus on:
1. Fixing the specific error mentioned above
2. Ensuring proper variable names and function calls
3. Using correct syntax and library functions
"""
    
    @staticmethod
    def get_fix_system_message(language: str) -> str:
        if language == 'zh':
            return "你是一个代码调试专家。"
        else:
            return "You are a code debugging expert."
    
    @staticmethod
    def get_solution_prompt(language: str, observation_text: str, 
                           has_computation_success: bool, has_visualization: bool) -> str:
        
        if language == 'zh':
            return f"""
以下是代码执行的结果：

<observation>
{observation_text}
</observation>

状态摘要：
- 计算：{'成功' if has_computation_success else '失败'}
- 可视化：{'已生成' if has_visualization else '未生成或失败'}

请在 <solution></solution> 标签内提供清晰的科学总结。

如果计算成功，请包括：
1. 对用户问题的直接回答
2. 关键数值结果及其解释
3. 如果生成了可视化，描述它显示的内容
4. 观察到的任何重要见解或模式

如果计算失败，请解释：
1. 出了什么问题
2. Agent为解决问题所做的努力，包括数据选择、工具选择和可视化尝试
3. 如何重新表述查询的建议
"""
        
        else:  # English
            return f"""
Here are the results of code execution:

<observation>
{observation_text}
</observation>

Status Summary:
- Computation: {'Success' if has_computation_success else 'Failed'}
- Visualization: {'Generated' if has_visualization else 'Not generated or failed'}

Please provide a clear scientific summary inside <solution></solution> tags.

If computation succeeded, include:
1. Direct answer to the user's question
2. Key numerical results and their interpretation
3. If visualization was generated, describe what it shows
4. Any important insights or patterns observed

If computation failed, explain:
1. What went wrong
2. Efforts made by the agent to solve the question, including data selection, tool choice, and visualization attempts
3. Suggestions for how to reformulate the query
"""
    
    @staticmethod
    def get_solution_input(language: str) -> str:
        if language == 'zh':
            return "生成最终中文回复。"
        else:
            return "Provide the final explanation in English."