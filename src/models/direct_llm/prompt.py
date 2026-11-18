"""
Prompts of the direct-LLM baseline: the question, the data and the output format.
"""
import json
from typing import Any, Dict

from src.models.prompt import PromptTemplates


class DirectLLMPrompts(PromptTemplates):
    """Code generation, repair and answer for DirectLLMAgent."""

    @staticmethod
    def get_code_generation_prompt(language: str, registry: Dict[str, Any], conventions: Dict[str, str],
                                   context_summary: str = "") -> str:
        """registry is every dataset's raw entry in configs/datasets.yaml; conventions the raw-file notes."""
        context_block = f"\n{context_summary}\n" if context_summary else ""
        if language == 'zh':
            return f"""你是一个数据分析助手。编写 Python 代码，用下面的数据集回答用户的问题。
{context_block}
本机可用的数据集，即每个数据集的登记条目（根目录、文件组织与命名、网格、时间轴、变量、各模式的坐标名与
单位，以及 notes）。部分 notes 描述的是某个 opener 或 loader 会做的处理；这里二者都不可用，notes 里说
opener 会做的事，需要你的代码自己完成：

{json.dumps(registry, indent=2, ensure_ascii=False)}

原始文件的跨数据集约定：

{json.dumps(conventions, indent=2, ensure_ascii=False)}

代码运行时已导入 numpy as np、pandas as pd、xarray as xr、matplotlib.pyplot as plt。

按以下两个代码块输出：

<compute>
# 计算：读取数据并计算结果。每个结果都存进变量；
# 答案依据的就是这个代码块结束时留下的变量。
</compute>

<visualize>
# 可选的绘图代码，可以使用 <compute> 中的变量。不要调用 plt.savefig()，
# 代码结束时仍打开的图会被保存。
</visualize>

惰性 (dask) 数组必须先 .compute()（或 float()）再存进要汇报的变量，否则只能看到占位符。
"""
        return f"""You are a data analysis assistant. Write Python code that answers the user's question with the
datasets below.
{context_block}
Datasets available on this machine: the registry entry of each (root directory, file layout and naming,
grid, time axis, variables, each model's coordinate names and units, and notes). Some notes describe
what an opener or loader does; neither is available here, so whatever a note says the opener does,
your code has to do itself:

{json.dumps(registry, indent=2, ensure_ascii=False)}

Cross-dataset conventions of the raw files:

{json.dumps(conventions, indent=2, ensure_ascii=False)}

The code runs with numpy as np, pandas as pd, xarray as xr and matplotlib.pyplot as plt already
imported.

Write it as two blocks:

<compute>
# Computation: read the data and compute the results. Store every result in a variable;
# the variables this block leaves behind are what the answer is written from.
</compute>

<visualize>
# Optional plotting code, which can use the variables from <compute>. Do not call
# plt.savefig(); the figures left open when it ends are saved.
</visualize>

Lazy (dask-backed) arrays must be computed (.compute() or float()) into the variables you want
reported, otherwise only a placeholder is seen.
"""

    @staticmethod
    def get_code_generation_input(language: str, user_query: str = "") -> str:
        if language == 'zh':
            return f"用户问题：{user_query}\n\n现在输出两个代码块。"
        return f"User question: {user_query}\n\nWrite the two code blocks now."

    @staticmethod
    def get_fix_prompt(language: str, code: str, exec_result: Dict, user_query: str,
                       full_specs: Dict, code_type: str,
                       upstream_code: str = "", available_variables: Dict = None) -> str:
        """
        full_specs is the registry, as in code generation. The code that ran before a
        visualization block and the variables it left are shown because the block
        runs in that namespace; without them its names cannot be repaired at all.
        """
        zh = language == 'zh'
        upstream_block = ""
        if upstream_code:
            upstream_block = (f"\n已经执行过的上游代码（它定义的变量在命名空间中可用）：\n<upstream_code>\n"
                              f"{upstream_code}\n</upstream_code>\n" if zh else
                              f"\nThe code that already ran before this block (its variables are in the "
                              f"namespace):\n<upstream_code>\n{upstream_code}\n</upstream_code>\n")
        variables_block = ""
        if available_variables:
            listing = json.dumps(available_variables, indent=2, ensure_ascii=False, default=str)
            variables_block = (f"\n命名空间中已有的变量（名称 -> 内容摘要）：\n{listing}\n" if zh else
                               f"\nVariables in the namespace (name -> summary of its contents):\n{listing}\n")
        registry = json.dumps(full_specs, indent=2, ensure_ascii=False)

        if zh:
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
{upstream_block}{variables_block}
本机可用的数据集：
{registry}

修复这个 {code_type} 代码块，并在 <execute></execute> 标签内输出修正后的完整版本。
"""
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
{upstream_block}{variables_block}
Datasets available on this machine:
{registry}

Fix this {code_type} block and output the corrected, complete version inside <execute></execute> tags.
"""

    @staticmethod
    def get_solution_prompt(language: str, observation_text: str,
                            has_computation_success: bool, has_visualization: bool,
                            user_query: str = "", plan: str = "", executed_code: str = "") -> str:
        """The shared answer step's signature; there is no plan in this baseline, so plan is unused."""
        if language == 'zh':
            return f"""
用户的问题：
{user_query}

执行的代码：
<executed_code>
{executed_code}
</executed_code>

状态：计算{'成功' if has_computation_success else '失败'}；可视化{'已生成' if has_visualization else '未生成或失败'}

结果：
<observation>
{observation_text}
</observation>

根据这些结果，在 <solution></solution> 标签内回答用户的问题。
"""
        return f"""
The user's question:
{user_query}

The code that ran:
<executed_code>
{executed_code}
</executed_code>

Status: computation {'succeeded' if has_computation_success else 'failed'}; visualization \
{'generated' if has_visualization else 'not generated or failed'}

Results:
<observation>
{observation_text}
</observation>

Answer the user's question inside <solution></solution> tags, based on these results.
"""
