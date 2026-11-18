"""
Prompts of the tool ablation (-T): the full system's prompts with the tools taken out.

The knowledge prompts (selection and condensing) are not overridden: retrieval
runs exactly as in the full system.
"""
from typing import Dict, List, Tuple

from src.models.prompt import PromptTemplates


Edit = Tuple[str, str]


def apply_edits(text: str, edits: List[Edit], prompt: str) -> str:
    """Apply each (old, new) edit once; an old text that is missing is an error."""
    for old, new in edits:
        if old not in text:
            raise ValueError(
                f"The tool-ablation edit for the {prompt} prompt no longer matches "
                f"src/models/prompt.py; update src/models/ablation_tools/prompt.py. "
                f"Missing text:\n{old}")
        text = text.replace(old, new, 1)
    return text


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

PLANNING_EDITS_EN: List[Edit] = [
    # The agent is built with no tools, so the overview renders as {}.
    ("Available tools:\n{}\n",
     "Available tools: none. There are no project functions in this setting; the analysis code opens\n"
     "the datasets' raw files itself with numpy, pandas and xarray.\n"),
    ("4. Tool selection: any question that needs data must include one of load_regular_data,\n"
     "   load_climatology_data or load_data_with_geolocator (the last one for named ocean regions); there\n"
     "   is no other way in. A published index series such as cpc_indices also goes through\n"
     "   load_regular_data, with 'product' choosing the series and 'variable' the region.\n",
     "4. Reading data: there are no loading functions. The code step opens the raw files of the datasets\n"
     "   you select and is shown where they live and how they are laid out; plan the reading as part of\n"
     "   the steps.\n"),
    ("5. Do not plan to compute area or volume weights by hand: select compute_basic_stat and use its\n"
     "   weighted operators.\n",
     "5. Regional and global means are area- or volume-weighted, never plain averages over lat/lon: plan\n"
     "   where the weights come from (for CMIP6, the model's own areacello / volcello fixed fields).\n"),
    ("Analyze the user's question and select appropriate data and tools.",
     "Analyze the user's question and select appropriate data."),
    ("<selected_tools>comma-separated list of tool names you will use</selected_tools>\n", ""),
    # The plan rules, shared by both knowledge arms.
    ("each naming the datasets, variables, models and tools it needs.",
     "each naming the datasets, variables and models it needs."),
    ("(the dataset and tool descriptions, and any domain knowledge)",
     "(the dataset descriptions, and any domain knowledge)"),
    ("3. Prefer the tools: when a listed tool performs a step, plan that tool and its operator rather than a\n"
     "   hand-written computation; plan your own code only for what no tool covers.\n",
     "3. There are no project tools to call: plan every step as an operation the code performs on the raw\n"
     "   data with the standard libraries (numpy, pandas, xarray).\n"),
    ("the datasets, variables, models and tools each step needs,",
     "the datasets, variables and models each step needs,"),
]

PLANNING_EDITS_EN_KNOWLEDGE: List[Edit] = [
    ("the\n   dataset and tool descriptions second,", "the\n   dataset descriptions second,"),
    ("(the formula, threshold, unit conversion, tool argument or check)",
     "(the formula, threshold, unit conversion or check)"),
    ("   severity should or info that end up unused may be left out.\n",
     "   severity should or info that end up unused may be left out.\n"
     "5. The notes may name functions of this project in backticks (`compute_basic_stat`,\n"
     "   `load_regular_data`, ...). None of them exists in this setting: read such a name as the operation\n"
     "   it stands for, and plan that operation on the raw data.\n"),
]

PLANNING_EDITS_ZH: List[Edit] = [
    ("可用工具：\n{}\n",
     "可用工具：无。本设置下没有项目函数，分析代码直接用 numpy、pandas、xarray 打开数据集的原始文件。\n"),
    ("4. 工具选择：凡是需要读数据，必须包含 load_regular_data / load_climatology_data /\n"
     "   load_data_with_geolocator 三者之一（命名海域用 load_data_with_geolocator），没有别的读取入口。\n"
     "   cpc_indices 这类已发布的指数序列也走 load_regular_data，用 product 选序列、variable 选区域。\n",
     "4. 读取数据：没有加载函数。代码步骤会直接打开所选数据集的原始文件，并会看到文件所在位置和组织方式；\n"
     "   请把读取作为计划的一部分。\n"),
    ("5. 区域平均、全球平均不要自己手算权重，选 compute_basic_stat 并使用其加权算子。\n",
     "5. 区域平均、全球平均必须按面积或体积加权，不能直接对 lat/lon 求平均：计划中写明权重从哪里来\n"
     "   （CMIP6 用模式自带的 areacello / volcello 固定场）。\n"),
    ("分析用户的问题并选择合适的数据和工具。", "分析用户的问题并选择合适的数据。"),
    ("<selected_tools>您将使用的工具名称，用逗号分隔</selected_tools>\n", ""),
    # The plan rules, shared by both knowledge arms.
    ("每一步写明需要哪些数据、变量、模式和工具；", "每一步写明需要哪些数据、变量和模式；"),
    ("（数据集说明、工具说明，以及领域知识——如果有的话）", "（数据集说明，以及领域知识——如果有的话）"),
    ("3. 优先使用工具：上面列出的工具能完成的步骤，计划用该工具及其算子，不要手写计算；只有没有工具覆盖的\n"
     "   部分才计划自己写代码。\n",
     "3. 本设置下没有项目工具可调用：每一步都规划为代码用常用库（numpy、pandas、xarray）在原始数据上完成的\n"
     "   操作。\n"),
    ("每一步需要哪些数据、变量、模式和工具，", "每一步需要哪些数据、变量和模式，"),
]

PLANNING_EDITS_ZH_KNOWLEDGE: List[Edit] = [
    ("其次是数据集说明与工具说明，", "其次是数据集说明，"),
    ("（公式、阈值、单位换算、工具参数或检查方法）", "（公式、阈值、单位换算或检查方法）"),
    ("   的片段若最终没有用上，可以不写。\n",
     "   的片段若最终没有用上，可以不写。\n"
     "5. 注意事项里可能在反引号中出现本项目的函数名（如 `compute_basic_stat`、`load_regular_data`）。本设置下\n"
     "   这些函数都不存在：把它理解为它所代表的操作，并在原始数据上规划这一操作。\n"),
]


# ---------------------------------------------------------------------------
# Code generation. full_specs carries the raw registry entries of the selected
# datasets in place of tool specs (AgentWithoutTools._code_reference).
# ---------------------------------------------------------------------------

CODE_EDITS_EN: List[Edit] = [
    ("Here are the full specifications for the tools you have selected:\n",
     "There are no project tools in this setting. Here is how the selected datasets are stored on disk:\n"
     "the registry entry of each (root directory, file layout and naming, grid, time axis, variables,\n"
     "each model's coordinate names and units, and notes). Some notes describe what this project's opener\n"
     "or loader does; neither exists here, so whatever a note says the opener does, your code has to do\n"
     "itself:\n"),
    ("Cross-dataset conventions (every opener already enforces these; use them to judge what can be\n"
     "compared directly):\n",
     "Cross-dataset conventions of the raw files (nothing enforces them for you; use them to judge what\n"
     "must be aligned before two products can be compared):\n"),
    (" — the IDs of the datasets above. Pass an ID string\n"
     "   straight to a loader's dataset argument (e.g. load_climatology_data(dataset=\"woa23\", ...)); never\n"
     "   pass a file path. Take variable and model names from the dataset entries above.\n",
     " — the IDs of the datasets above. There are no\n"
     "   loaders: open the raw files yourself under each dataset's root, following the file layout in its\n"
     "   entry. Take variable and model names from the dataset entries above.\n"),
    ("# Call a loading tool to read the data, then process it according to the user's question\n",
     "# Open the raw data files, then process them according to the user's question\n"),
    ("- The loaders return LAZY (dask-backed) xr.DataArray objects. Every quantity you want reported must\n"
     "  be materialised explicitly with .compute() (or float()) into its own small variable, otherwise the\n"
     "  downstream step sees a dask placeholder and no numbers at all. Large intermediates can stay lazy.\n",
     "- Arrays opened with chunks (xr.open_dataset(..., chunks=...), xr.open_mfdataset) are LAZY\n"
     "  (dask-backed). Every quantity you want reported must be materialised explicitly with .compute() (or\n"
     "  float()) into its own small variable, otherwise the downstream step sees a dask placeholder and no\n"
     "  numbers at all. Large intermediates can stay lazy.\n"),
    ("- For regional or global means use compute_basic_stat with its weighted operators\n"
     "  (area_weighted_mean / volume_weighted_mean). Do not average over lat/lon with a plain .mean(), and\n"
     "  do not fetch areacello by hand.\n",
     "- Regional or global means are area- or volume-weighted: for a CMIP6 model by its own areacello /\n"
     "  volcello fixed fields, for a regular latitude-longitude grid by the area of each cell. Do not\n"
     "  average over lat/lon with a plain .mean().\n"),
    ("Subset them through the loader's lat_range /\n  lon_range arguments instead.\n",
     "Subset them with a mask built from their 2-D\n  latitude/longitude coordinates instead.\n"),
    ("convert with transform_unit before comparing them.\n",
     "convert with the in situ density (mol m-3 = umol kg-1 * rho / 1000) before\n  comparing them.\n"),
    ("- For a region that crosses the dateline (the equatorial Pacific, 120E-280E), just pass lon_range to\n"
     "  the loading tool: it returns a continuous longitude axis that the weighted statistics and\n"
     "  derived-quantity tools accept as it is. Do not assemble such a region by hand with %360, by\n"
     "  concatenating two sel windows, or by sortby — an array rebuilt that way has left the dataset's own\n"
     "  grid, and weighting it fails with \"these lon values have no cell on the dataset's own grid\".\n",
     "- For a region that crosses the dateline (the equatorial Pacific, 120E-280E), mind each product's\n"
     "  native longitude convention (see the conventions above): such a window is contiguous on a 0-360\n"
     "  axis but wraps around on a -180-180 one.\n"),
    ("- Arrays from the loaders carry provenance attributes (source_dataset, grid_type, ...) that the\n"
     "  weighted statistics and derived-quantity tools use to find the right grid. The execution\n"
     "  environment runs with keep_attrs enabled, so ordinary xarray operations preserve them; but when a\n"
     "  step changes the physical quantity (a ratio, a difference of unlike things), set the new units\n"
     "  yourself with .assign_attrs(units=...), or the reported units will be wrong.\n",
     "- The execution environment runs with keep_attrs enabled, so ordinary xarray operations preserve\n"
     "  attributes; but when a step changes the physical quantity (a ratio, a difference of unlike\n"
     "  things), set the new units yourself with .assign_attrs(units=...), or the reported units will be\n"
     "  wrong.\n"),
    ("- A field on a curvilinear grid has 2-D coordinates that are NOT called lat/lon. Read the names off\n"
     "  the array with `da.attrs.get(\"lat_coord\", \"lat\")` / `da.attrs.get(\"lon_coord\", \"lon\")` rather than\n"
     "  indexing `da.coords[\"lat\"]`, which raises KeyError on CESM2 and MPI-ESM1-2-LR.\n",
     "- A field on a curvilinear grid has 2-D coordinates whose names differ by model. Take them from the\n"
     "  model's lat_coord / lon_coord in the dataset entry above rather than assuming `lat` / `lon`.\n"),
    ("first put it on a\n  regular grid with `regrid_horizontal(field, target=\"woa23\")` (for display only;",
     "first put it on a\n  regular 1-degree longitude-latitude grid, interpolating from its 2-D latitude/longitude (e.g. with\n"
     "  scipy.interpolate) (for display only;"),
]

CODE_EDITS_ZH: List[Edit] = [
    ("以下是您选择的工具的完整规格说明：\n",
     "本设置下没有项目工具。以下是所选数据集在磁盘上的存储方式，即每个数据集的登记条目（根目录、文件组织与\n"
     "命名、网格、时间轴、变量、各模式的坐标名与单位，以及 notes）。部分 notes 描述的是本项目的 opener 或\n"
     "loader 会做的处理；这里二者都不存在，notes 里说 opener 会做的事，需要你的代码自己完成：\n"),
    ("跨数据集约定（所有 opener 已统一遵守，编码时据此判断能否直接比较）：\n",
     "原始文件的跨数据集约定（没有任何环节替你统一，编码时据此判断两个产品比较前需要先对齐什么）：\n"),
    ("，即上面这些数据集的 ID。\n"
     "   把 ID 字符串直接传给加载工具的 dataset 参数（如 load_climatology_data(dataset=\"woa23\", ...)），\n"
     "   不要使用文件路径。变量名与模式名按上面的数据集说明填写。\n",
     "，即上面这些数据集的 ID。\n"
     "   没有加载函数：按各数据集条目中的文件组织方式，自行在其根目录下打开原始文件。变量名与模式名按上面的\n"
     "   数据集说明填写。\n"),
    ("# 调用加载工具读取数据，之后根据用户问题处理数据\n",
     "# 打开原始数据文件，之后根据用户问题处理数据\n"),
    ("- 加载工具返回的是惰性 (dask) 的 xr.DataArray。",
     "- 带 chunks 打开的数组（xr.open_dataset(..., chunks=...)、xr.open_mfdataset）是惰性 (dask) 的。"),
    ("- 区域/全球平均用 compute_basic_stat 的加权算子（area_weighted_mean / volume_weighted_mean），\n"
     "  不要用 .mean() 直接对 lat/lon 求平均，也不要手工取 areacello。\n",
     "- 区域/全球平均要按面积或体积加权：CMIP6 模式用其自带的 areacello / volcello 固定场，规则经纬网格用\n"
     "  每个格点的面积。不要用 .mean() 直接对 lat/lon 求平均。\n"),
    ("按经纬度取子集要交给加载工具的 lat_range/lon_range 参数。\n",
     "按经纬度取子集要用其二维经纬度坐标构造掩膜。\n"),
    ("互相比较前先用\n  transform_unit 换算。\n",
     "互相比较前先用现场密度\n  换算（mol m-3 = umol kg-1 * rho / 1000）。\n"),
    ("- 跨日界线的区域（如赤道太平洋 120E-280E）直接把 lon_range 传给加载工具即可，它会返回连续的经度轴，\n"
     "  加权统计与派生量工具都能接着用。不要自己用 %360、sel 拼接或 sortby 去手工拼跨日界线的区域，那样得到\n"
     "  的数组已经脱离数据集原始网格，加权时会报「这些 lon 值不在数据集网格上」。\n",
     "- 跨日界线的区域（如赤道太平洋 120E-280E）要注意各产品原生的经度约定（见上面的约定）：这样的窗口在\n"
     "  0-360 经度轴上是连续的，在 -180-180 经度轴上则会跨过边界。\n"),
    ("- 加载工具返回的数组带有 source_dataset / grid_type 等溯源属性，加权统计与派生量工具靠它们找网格。\n"
     "  执行环境已开启 keep_attrs，普通 xarray 运算不会丢掉这些属性；",
     "- 执行环境已开启 keep_attrs，普通 xarray 运算不会丢掉属性；"),
    ("- 曲线网格上的场，其二维坐标并不叫 lat/lon。请用 `da.attrs.get(\"lat_coord\", \"lat\")` /\n"
     "  `da.attrs.get(\"lon_coord\", \"lon\")` 读出坐标名，不要直接索引 `da.coords[\"lat\"]`——\n"
     "  在 CESM2 和 MPI-ESM1-2-LR 上会抛 KeyError。\n",
     "- 曲线网格上的场，其二维坐标名因模式而异。请按上面数据集条目中该模式的 lat_coord / lon_coord 取坐标名，\n"
     "  不要默认它们叫 `lat` / `lon`。\n"),
    ("先用\n  `regrid_horizontal(field, target=\"woa23\")` 放到规则网格上（仅用于显示，",
     "先从其二维经纬度\n  插值到 1° 规则经纬网格上（例如用 scipy.interpolate；仅用于显示，"),
]


# ---------------------------------------------------------------------------
# Repair and final answer
# ---------------------------------------------------------------------------

FIX_EDITS_EN: List[Edit] = [
    ("Available tools and their specifications:\n",
     "How the selected datasets are stored on disk (their registry entries; there are no project tools):\n"),
]

FIX_EDITS_ZH: List[Edit] = [
    ("可用工具及其规格说明：\n", "所选数据集在磁盘上的存储方式（登记条目；没有项目工具）：\n"),
]

SOLUTION_EDITS_EN: List[Edit] = [
    ("must match the tools and arguments the code actually\ncalled.",
     "must match the functions and arguments the code actually\ncalled."),
    ("including data selection, tool choice, computation and visualization attempts",
     "including data selection, computation and visualization attempts"),
]

SOLUTION_EDITS_ZH: List[Edit] = [
    ("一律以代码里实际调用的工具和参数为准", "一律以代码里实际调用的函数和参数为准"),
    ("包括数据选择、工具选择以及计算与可视化尝试", "包括数据选择以及计算与可视化尝试"),
]


class PromptTemplatesWithoutTools(PromptTemplates):
    """The full system's prompts with every tool reference edited out (see the module docstring)."""

    @staticmethod
    def get_planning_prompt(language: str, tool_overview: Dict, data_overview: Dict, context_summary: str,
                            knowledge: str = "") -> str:
        text = PromptTemplates.get_planning_prompt(language, tool_overview, data_overview,
                                                   context_summary, knowledge=knowledge)
        if language == 'zh':
            edits = PLANNING_EDITS_ZH + (PLANNING_EDITS_ZH_KNOWLEDGE if knowledge else [])
        else:
            edits = PLANNING_EDITS_EN + (PLANNING_EDITS_EN_KNOWLEDGE if knowledge else [])
        return apply_edits(text, edits, "planning")

    @staticmethod
    def get_code_generation_prompt(language: str, full_specs: Dict, user_query: str,
                                   dataset_specs: Dict, conventions: Dict, plan: str = "") -> str:
        text = PromptTemplates.get_code_generation_prompt(language, full_specs, user_query,
                                                          dataset_specs, conventions, plan=plan)
        return apply_edits(text, CODE_EDITS_ZH if language == 'zh' else CODE_EDITS_EN, "code-generation")

    @staticmethod
    def get_fix_prompt(language: str, code: str, exec_result: Dict, user_query: str,
                       full_specs: Dict, code_type: str,
                       upstream_code: str = "", available_variables: Dict = None) -> str:
        text = PromptTemplates.get_fix_prompt(language, code, exec_result, user_query, full_specs,
                                              code_type, upstream_code=upstream_code,
                                              available_variables=available_variables)
        return apply_edits(text, FIX_EDITS_ZH if language == 'zh' else FIX_EDITS_EN, "fix")

    @staticmethod
    def get_solution_prompt(language: str, observation_text: str,
                            has_computation_success: bool, has_visualization: bool,
                            user_query: str = "", plan: str = "", executed_code: str = "") -> str:
        text = PromptTemplates.get_solution_prompt(language, observation_text, has_computation_success,
                                                   has_visualization, user_query=user_query, plan=plan,
                                                   executed_code=executed_code)
        return apply_edits(text, SOLUTION_EDITS_ZH if language == 'zh' else SOLUTION_EDITS_EN, "solution")
