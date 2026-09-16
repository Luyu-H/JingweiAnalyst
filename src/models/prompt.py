import json
from typing import Dict, Any


class PromptTemplates:
    """Collection of all prompt templates for the agent."""     
    
    @staticmethod
    def get_planning_prompt(language: str, tool_overview: Dict, data_overview: Dict, context_summary: str) -> str:

        if language == 'zh':
            return f"""你是一个海洋科学数据分析助手，支持对话历史功能。

{context_summary}

可用工具：
{json.dumps(tool_overview, indent=2, ensure_ascii=False)}

可用数据集：
{json.dumps(data_overview, indent=2, ensure_ascii=False)}

数据选择规则：
1. <selected_data> 里只能填上面出现的数据集 ID（如 woa23、cmip6_historical），不要填文件路径，
   也不要自己拼接变量名。一个数据集 ID 下含多个变量，CMIP6 还含多个模式，这些在下一步选择。
2. 问长期趋势、年际变率、情景预估 -> 选 timeseries 类数据集（cmip6_historical、cmip6_ssp245、ersst）；
   问平均态、季节循环、垂向结构 -> 选 climatology 类（woa23、glodap、seanoe_mld）；
   问 ENSO 监测指数本身 -> 选 cpc_indices。
3. 模式与观测对比、多变量联合分析时，可以同时选多个数据集，用逗号分隔。
4. 工具选择：凡是需要读数据，必须包含 load_regular_data / load_climatology_data /
   load_data_with_geolocator 三者之一（命名海域用 load_data_with_geolocator）。open_* 系列是
   底层读取器，只在没有对应 load_* 入口时才直接用（例如 cpc_indices 只能用 open_cpc_index）。
5. 区域平均、全球平均不要自己手算权重，选 compute_basic_stat 并使用其加权算子。

分析用户的问题并选择合适的数据和工具。严格按照以下格式输出：
<plan>关于需要哪些数据、变量、模式和工具的逐步推理</plan>
<selected_data>您将使用的数据集 ID，用逗号分隔</selected_data>
<selected_tools>您将使用的工具名称，用逗号分隔</selected_tools>

尤其需要注意，回复中必须包含所需的所有XML标签，并且标签必须正确闭合。
"""
        
        else:  # English
            return f"""You are an oceanographic data analysis assistant with conversation history support.

{context_summary}

Available tools:
{json.dumps(tool_overview, indent=2)}

Available datasets:
{json.dumps(data_overview, indent=2)}

Data selection rules:
1. <selected_data> must contain dataset IDs exactly as listed above (e.g. woa23, cmip6_historical).
   Never put a file path there, and never invent a name. One dataset ID holds many variables, and
   CMIP6 also holds several models; those are chosen in the next step, not here.
2. Long-term trends, interannual variability and scenario projections need a timeseries dataset
   (cmip6_historical, cmip6_ssp245, ersst). Mean state, seasonal cycle and vertical structure need a
   climatology (woa23, glodap, seanoe_mld). A question about the ENSO monitoring index itself needs
   cpc_indices.
3. Select several datasets, comma-separated, when the question compares a model against observations
   or needs variables that live in different products.
4. Tool selection: any question that needs data must include one of load_regular_data,
   load_climatology_data or load_data_with_geolocator (the last one for named ocean regions). The
   open_* tools are the low-level readers; use one directly only when no load_* entry point covers
   the dataset (cpc_indices is reachable only through open_cpc_index).
5. Do not plan to compute area or volume weights by hand: select compute_basic_stat and use its
   weighted operators.

Analyze the user's question and select appropriate data and tools. Follow this format strictly:
<plan>step-by-step reasoning about which datasets, variables, models and tools are needed</plan>
<selected_data>comma-separated list of dataset IDs you will use</selected_data>
<selected_tools>comma-separated list of tool names you will use</selected_tools>

IMPORTANT: Your response MUST contain all required XML tags, and tags must be properly closed.
"""
    
    @staticmethod
    def get_planning_input(language: str, user_query: str) -> str:
        if language == 'zh':
            return f"用户问题：{user_query}\n\n现在输出您的推理和选择。"
        else:
            return f"User question: {user_query}\n\nNow output your reasoning and selections."
    
    @staticmethod
    def get_code_generation_prompt(language: str, full_specs: Dict, user_query: str,
                                   dataset_specs: Dict, conventions: Dict, plan: str = "") -> str:
        """
        Build the code-generation system prompt.

        dataset_specs carries the full registry entry of every dataset the planning
        step picked (variables, models, time span, climatology period), and
        conventions carries the registry's cross-dataset rules. Without them this
        step would only see tool signatures and would have to guess which model or
        variable names actually exist.
        """
        dataset_ids = list(dataset_specs.keys())

        if language == 'zh':
            plan_block = f"""
上一步的分析计划：
{plan}
""" if plan else ""

            return f"""以下是您选择的工具的完整规格说明：

{json.dumps(full_specs, indent=2, ensure_ascii=False)}

以下是您选择的数据集的完整说明（变量、模式、时间范围、气候基准期）：

{json.dumps(dataset_specs, indent=2, ensure_ascii=False)}

跨数据集约定（所有 opener 已统一遵守，编码时据此判断能否直接比较）：

{json.dumps(conventions, indent=2, ensure_ascii=False)}
{plan_block}
生成两个独立的代码块来解决："{user_query}"

环境变量与可用库：
1. `selected_datasets` = {json.dumps(dataset_ids, ensure_ascii=False)}，即上面这些数据集的 ID。
   把 ID 字符串直接传给加载工具的 dataset 参数（如 load_climatology_data(dataset="woa23", ...)），
   不要使用文件路径。变量名与模式名按上面的数据集说明填写。
2. 可用库（已注入，无需 import）：numpy as np, pandas as pd, xarray as xr, matplotlib.pyplot as plt

按以下格式生成代码：

<compute>
# 仅计算代码
# 调用加载工具读取数据，之后根据用户问题处理数据
# 执行计算并将结果存储在变量中，比如 result, mean_value, trend 等
</compute>

<visualize>
# 仅可视化代码（此部分可选）
# 使用计算结果（它们将作为变量可用）
# 创建带有适当标签和标题的matplotlib图表
# 使用 plt.figure(), plt.plot(), plt.xlabel() 等
# 不要使用 plt.savefig()
</visualize>

重要提示：
- 计算代码必须完整且独立
- 加载工具返回的是惰性 (dask) 的 xr.DataArray。最终要汇报的结论量必须显式 .compute()（或 float()）
  求值后存进单独的小变量，否则下游只能看到 dask 占位符而读不到数值。中间大数组不必 compute。
- 区域/全球平均用 compute_basic_stat 的加权算子（area_weighted_mean / volume_weighted_mean），
  不要用 .mean() 直接对 lat/lon 求平均，也不要手工取 areacello。
- CESM2 与 MPI-ESM1-2-LR 是曲线网格，其索引维 (nlat/nlon, j/i) 不是经纬度，不能用 .sel() 按值切片；
  按经纬度取子集要交给加载工具的 lat_range/lon_range 参数。
- CMIP6 浓度是 mol m-3（体积基），WOA23 与 GLODAP 是 umol kg-1（质量基），互相比较前先用
  transform_unit 换算。
- 跨日界线的区域（如赤道太平洋 120E-280E）直接把 lon_range 传给加载工具即可，它会返回连续的经度轴，
  加权统计与派生量工具都能接着用。不要自己用 %360、sel 拼接或 sortby 去手工拼跨日界线的区域，那样得到
  的数组已经脱离数据集原始网格，加权时会报「这些 lon 值不在数据集网格上」。
- 加载工具返回的数组带有 source_dataset / grid_type 等溯源属性，加权统计与派生量工具靠它们找网格。
  执行环境已开启 keep_attrs，普通 xarray 运算不会丢掉这些属性；但若某步运算改变了物理量纲
  （如求比值、做差分），请自行用 .assign_attrs(units=...) 更新 units，否则单位标注会失真。
- 可视化代码可以引用计算中的变量
- 只要结果还有可展示的维度就必须给出 <visualize>：时间序列、季节循环、垂向廓线、空间分布图、
  多组对比，都属于此列，即使用户没有明说要图也要画。只有当答案就是一个标量或一张几行的表格、
  确实无图可画时，才把 <visualize> 留空。
"""

        else:  # English
            plan_block = f"""
Plan from the previous step:
{plan}
""" if plan else ""

            return f"""Here are the full specifications for the tools you have selected:

{json.dumps(full_specs, indent=2)}

Here are the full entries for the datasets selected in the previous step (variables, models, time
span, climatology period):

{json.dumps(dataset_specs, indent=2, ensure_ascii=False)}

Cross-dataset conventions (every opener already enforces these; use them to judge what can be
compared directly):

{json.dumps(conventions, indent=2, ensure_ascii=False)}
{plan_block}
Generate TWO separate code blocks to solve: "{user_query}"

Context variables and available libraries:
1. `selected_datasets` = {json.dumps(dataset_ids)} — the IDs of the datasets above. Pass an ID string
   straight to a loader's dataset argument (e.g. load_climatology_data(dataset="woa23", ...)); never
   pass a file path. Take variable and model names from the dataset entries above.
2. Available libraries (already injected, no import needed): numpy as np, pandas as pd,
   xarray as xr, matplotlib.pyplot as plt

Generate code in this format:

<compute>
# Computation code only
# Call a loading tool to read the data, then process it according to the user's question
# Perform calculations and store results in variables (e.g. result, mean_value, trend, etc.)
</compute>

<visualize>
# Visualization code only (this section is optional)
# Use results from computation (they will be available as variables)
# Create matplotlib plots with proper labels and titles
# Use plt.figure(), plt.plot(), plt.xlabel(), etc.
# DO NOT use plt.savefig()
</visualize>

Important:
- Computation code MUST be complete and standalone
- The loaders return LAZY (dask-backed) xr.DataArray objects. Every quantity you want reported must
  be materialised explicitly with .compute() (or float()) into its own small variable, otherwise the
  downstream step sees a dask placeholder and no numbers at all. Large intermediates can stay lazy.
- For regional or global means use compute_basic_stat with its weighted operators
  (area_weighted_mean / volume_weighted_mean). Do not average over lat/lon with a plain .mean(), and
  do not fetch areacello by hand.
- CESM2 and MPI-ESM1-2-LR are on curvilinear grids: their index dimensions (nlat/nlon, j/i) are NOT
  geographic and cannot be sliced by value with .sel(). Subset them through the loader's lat_range /
  lon_range arguments instead.
- CMIP6 concentrations are in mol m-3 (volume basis) while WOA23 and GLODAP are in umol kg-1 (mass
  basis); convert with transform_unit before comparing them.
- For a region that crosses the dateline (the equatorial Pacific, 120E-280E), just pass lon_range to
  the loading tool: it returns a continuous longitude axis that the weighted statistics and
  derived-quantity tools accept as it is. Do not assemble such a region by hand with %360, by
  concatenating two sel windows, or by sortby — an array rebuilt that way has left the dataset's own
  grid, and weighting it fails with "these lon values have no cell on the dataset's own grid".
- Arrays from the loaders carry provenance attributes (source_dataset, grid_type, ...) that the
  weighted statistics and derived-quantity tools use to find the right grid. The execution
  environment runs with keep_attrs enabled, so ordinary xarray operations preserve them; but when a
  step changes the physical quantity (a ratio, a difference of unlike things), set the new units
  yourself with .assign_attrs(units=...), or the reported units will be wrong.
- Visualization code can reference variables from computation
- Provide a <visualize> block whenever the result still has something to show — a time series, a
  seasonal cycle, a depth profile, a spatial map, or several groups to compare — even when the
  question does not ask for a figure in so many words. Leave it empty only when the answer is a
  single number or a short table with nothing to plot.
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
                           has_computation_success: bool, has_visualization: bool,
                           user_query: str = "", plan: str = "", executed_code: str = "") -> str:
        """
        Build the final-answer prompt.

        The question, the plan and the code that actually ran are all included, so
        that the methodology is reported from the code. Inferred from the names of
        the result variables alone, a baseline, a weighting or a unit is easy to
        state wrongly, and a wrong methodology invalidates the whole summary.
        """
        if language == 'zh':
            return f"""
用户的问题：
{user_query}

分析计划：
{plan}

实际执行并产生下列结果的代码：
<executed_code>
{executed_code}
</executed_code>

以下是代码执行的结果：

状态摘要：
- 计算：{'成功' if has_computation_success else '失败'}
- 可视化：{'已生成' if has_visualization else '未生成或失败'}

<observation>
{observation_text}
</observation>

请在 <solution></solution> 标签内提供清晰的科学总结。

方法学描述必须严格依据上面的 <executed_code>，不要依据变量名猜测：基准期、气候态来源、加权方式、
单位、区域与时间范围、滞后范围等，一律以代码里实际调用的工具和参数为准。代码里没做过的事不要写。
观测结果中的变量名可能有歧义，遇到歧义时以代码为准。

如果计算成功，请包括：
1. 对用户问题的直接回答
2. 关键数值结果及其解释
3. 如果生成了可视化，描述它显示的内容。注意描述可视化逻辑和方法，而不是解释具体的视觉结果。
4. 观察到的任何重要见解或模式

如果计算失败，请解释：
1. 出了什么问题
2. Agent为解决问题所做的努力，包括数据选择、工具选择以及计算与可视化尝试
3. 如何重新表述查询的建议
"""
        
        else:  # English
            return f"""
The user's question:
{user_query}

The analysis plan:
{plan}

The code that actually ran and produced the results below:
<executed_code>
{executed_code}
</executed_code>

Here are the results of code execution:

Status Summary:
- Computation: {'Success' if has_computation_success else 'Failed'}
- Visualization: {'Generated' if has_visualization else 'Not generated or failed'}

<observation>
{observation_text}
</observation>

Please provide a clear scientific summary inside <solution></solution> tags.

Describe the methodology strictly from the <executed_code> above, never by inferring it from variable
names: the baseline period, the source of any climatology, the weighting, the units, the region and
time range, the range of lags — all of it must match the tools and arguments the code actually
called. Do not report a step the code did not take. Where a variable name in the observation is
ambiguous, the code decides what it means.

If computation succeeded, include:
1. Direct answer to the user's question
2. Key numerical results and their interpretation
3. If visualization was generated, describe what it shows. The description should outline the visualization's intended logic and methodology, without interpreting the specific visual result.
4. Any important insights or patterns observed

If computation failed, explain:
1. What went wrong
2. Efforts made by the agent to solve the question, including data selection, tool choice, computation and visualization attempts
3. Suggestions for how to reformulate the query
"""
    
    @staticmethod
    def get_solution_input(language: str) -> str:
        if language == 'zh':
            return "生成最终中文回复。"
        else:
            return "Provide the final explanation in English."