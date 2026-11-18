import json
from typing import Dict


class PromptTemplates:
    """Collection of all prompt templates for the agent."""     
    
    @staticmethod
    def get_planning_prompt(language: str, tool_overview: Dict, data_overview: Dict, context_summary: str,
                            knowledge: str = "") -> str:
        """
        Build the planning system prompt.
        """
        if language == 'zh':
            knowledge_block = f"""
为本问题检索到的领域知识（按 [片段id] (severity) 分组）：
<domain_knowledge>
{knowledge}
</domain_knowledge>
""" if knowledge else ""
            plan_rules = """
计划的写法：
1. 计划写成编号步骤，每一步写明需要哪些数据、变量、模式和工具；最后写一节 "Other notes:"，每条一行，
   只放答案必须声明的注意事项。
2. 不要自行添加计算约束。步骤里用到的每个定义、权重、掩膜、阈值、插值、筛选或订正，都必须来自问题，
   或来自上面给出的材料（数据集说明、工具说明，以及领域知识——如果有的话）；它们都没有写的，不要自己
   设计（例如额外的掩膜、按部分层厚加权、筛选阈值）。
3. 优先使用工具：上面列出的工具能完成的步骤，计划用该工具及其算子，不要手写计算；只有没有工具覆盖的
   部分才计划自己写代码。
4. 中间量不要四舍五入：后面的步骤使用前面步骤未经取整的结果，只在报告数值时按问题要求的精度取整。
5. 步骤、步骤下的注意事项和 Other notes 全部写在 <plan> 标签内；标签外的内容不会传给后续步骤。
"""
            knowledge_rules = """
领域知识的使用规则：
1. 优先级：用户在问题里明确给出的定义和步骤最优先，其次是数据集说明与工具说明，最后才是知识片段。
   知识片段用来决定问题没有规定的选择、找出需要说明的注意事项；某个片段与问题规定的做法冲突时，
   照问题的做法执行，并把差异写进 Other notes，让最终答案能够解释。
2. 在某条知识【第一次起作用】的那一步下面加一行 "- 注意 [片段id]：……"，写出它能直接
   执行的内容（公式、阈值、单位换算、工具参数或检查方法）。整个计划中每个片段只能有【一行】这样的注意
   事项（在某一步下或在 Other notes 中，二者选一）。后面的步骤若也依赖它，只在该步的正文里写上
   "[片段id]"，不要再加注意事项行，也不要复述内容——下一步会读到整个计划。不要把知识原文抄进计划。
3. Other notes 中还要写出问题定义与协议定义的差异；来自知识片段的条目格式为 "- [片段id] ……"。
4. severity 为 must 的片段都必须至少引用一次（在某一步或 Other notes 中）。severity 为 should 或 info
   的片段若最终没有用上，可以不写。
""" if knowledge else ""
            plan_format = ("<plan>编号步骤：每一步需要哪些数据、变量、模式和工具，步骤下列出对应的知识注意事项，"
                           "最后是 Other notes</plan>" if knowledge else
                           "<plan>编号步骤：每一步需要哪些数据、变量、模式和工具，最后是 Other notes</plan>")

            return f"""你是一个海洋科学数据分析助手，支持对话历史功能。

{context_summary}

可用工具：
{json.dumps(tool_overview, indent=2, ensure_ascii=False)}

可用数据集：
{json.dumps(data_overview, indent=2, ensure_ascii=False)}
{knowledge_block}
数据选择规则：
1. <selected_data> 里只能填上面出现的数据集 ID（如 woa23、cmip6_historical），不要填文件路径，
   也不要自己拼接变量名。一个数据集 ID 下含多个变量，CMIP6 还含多个模式，这些在下一步选择。
2. 问长期趋势、年际变率、情景预估 -> 选 timeseries 类数据集（cmip6_historical、cmip6_ssp245、ersst）；
   问平均态、季节循环、垂向结构 -> 选 climatology 类（woa23、glodap、seanoe_mld）；
   问 ENSO 监测指数本身 -> 选 cpc_indices。
3. 模式与观测对比、多变量联合分析时，可以同时选多个数据集，用逗号分隔。
4. 工具选择：凡是需要读数据，必须包含 load_regular_data / load_climatology_data /
   load_data_with_geolocator 三者之一（命名海域用 load_data_with_geolocator），没有别的读取入口。
   cpc_indices 这类已发布的指数序列也走 load_regular_data，用 product 选序列、variable 选区域。
5. 区域平均、全球平均不要自己手算权重，选 compute_basic_stat 并使用其加权算子。
{plan_rules}{knowledge_rules}
分析用户的问题并选择合适的数据和工具。严格按照以下格式输出：
{plan_format}
<selected_data>您将使用的数据集 ID，用逗号分隔</selected_data>
<selected_tools>您将使用的工具名称，用逗号分隔</selected_tools>

尤其需要注意，回复中必须包含所需的所有XML标签，并且标签必须正确闭合。
"""
        
        else:  # English
            knowledge_block = f"""
Domain knowledge retrieved for this question (grouped by [fragment-id] (severity)):
<domain_knowledge>
{knowledge}
</domain_knowledge>
""" if knowledge else ""
            plan_rules = """
How to write the plan:
1. Write the plan as numbered steps, each naming the datasets, variables, models and tools it needs.
   End it with an "Other notes:" section, one line each, holding only the caveats the answer must
   state.
2. Add no computational constraint of your own. Every definition, weighting, mask, threshold,
   interpolation, filter or correction a step applies must come from the question or from what is
   given above (the dataset and tool descriptions, and any domain knowledge); do not devise one that
   none of them states (an extra mask, a partial-layer weighting, a screening threshold).
3. Prefer the tools: when a listed tool performs a step, plan that tool and its operator rather than a
   hand-written computation; plan your own code only for what no tool covers.
4. Do not round intermediate values: later steps use the unrounded results of earlier ones, and a
   value is rounded only where it is reported, to the precision the question asks for.
5. The steps, any notes under them and Other notes all go inside the <plan> tags; nothing outside them
   is passed on to the later steps.
"""
            knowledge_rules = """
How to use the domain knowledge:
1. Priority: definitions and procedures the user states explicitly in the question come first, the
   dataset and tool descriptions second, the knowledge fragments last. Use the fragments to settle
   choices the question leaves open and to find the caveats that must be stated. When a fragment
   conflicts with what the question prescribes, follow the question and record the difference under
   Other notes so the answer can explain it.
2. Under the step where a piece of knowledge FIRST matters, add one line "- Note [fragment-id]: ..."
   carrying what to do (the formula, threshold, unit conversion, tool argument or check). Each
   fragment gets exactly ONE such note line in the whole plan, under a step or under Other notes, not
   both. A later step that also depends on it names "[fragment-id]" inline in its own text, with no
   further note line and no restatement - the next stage reads the whole plan. Do not copy the
   knowledge text into the plan.
3. Other notes also hold the differences between the question's definition and the protocol; a line
   that comes from a fragment takes the form "- [fragment-id] ...".
4. Cite every fragment of severity must at least once, in a step or under Other notes. Fragments of
   severity should or info that end up unused may be left out.
""" if knowledge else ""
            plan_format = ("<plan>numbered steps: the datasets, variables, models and tools each step needs, "
                           "with the knowledge notes under each step, then Other notes</plan>" if knowledge else
                           "<plan>numbered steps: the datasets, variables, models and tools each step needs, "
                           "then Other notes</plan>")

            return f"""You are an oceanographic data analysis assistant with conversation history support.

{context_summary}

Available tools:
{json.dumps(tool_overview, indent=2)}

Available datasets:
{json.dumps(data_overview, indent=2)}
{knowledge_block}
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
   load_climatology_data or load_data_with_geolocator (the last one for named ocean regions); there
   is no other way in. A published index series such as cpc_indices also goes through
   load_regular_data, with 'product' choosing the series and 'variable' the region.
5. Do not plan to compute area or volume weights by hand: select compute_basic_stat and use its
   weighted operators.
{plan_rules}{knowledge_rules}
Analyze the user's question and select appropriate data and tools. Follow this format strictly:
{plan_format}
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
    def get_knowledge_selection_prompt(language: str, candidates, context_summary: str,
                                       max_selected: int) -> str:
        """
        Build the system prompt for choosing knowledge fragments before planning.
        """
        rows = []
        for c in candidates:
            matched = "; ".join(f"{k}: {', '.join(map(str, v))}" for k, v in c.get("signals", {}).items())
            rows.append(f"- {c['id']} | {c['type']} | {c['severity']} | score {c['score']}"
                        f"{' | matched ' + matched if matched else ''}\n  {c['summary']}")
        catalogue = "\n".join(rows)
        context_block = f"\n{context_summary}\n" if context_summary else ""

        if language == 'zh':
            return f"""你负责在制定分析计划之前，从项目知识库中挑出与用户问题相关的知识片段。
{context_block}
候选片段（id | 类型 | severity | 词法预筛分数 | 命中的词；下一行是摘要）：
{catalogue}

severity 含义：must = 忽略会算出错误的数值；should = 忽略会导致误导性的比较或解读；info = 解读结果的背景。

选择规则：
1. 选可能改变这道题的做法或结论表述的片段：它约束某个计算步骤（单位、加权、定义、网格、时间平均的
   顺序），或者答案必须声明它带来的注意事项。
2. 依据摘要判断片段是否适用，不要依据词语重合。预筛分数只是参考：问题提到某个词，不代表那个片段适用。
3. 这一步的目标是不漏：下一步会通读所选片段的全文，删掉不适用的部分。拿不准某个 must 或 should 片段
   是否适用时，选上它；所有数据处理都会遇到的通用事项，排在针对本题变量、数据集和方法的片段之后。
4. 最多选 {max_selected} 个。一个都不适用时，输出空标签。
5. 不准编造：只能从上面的候选列表中选 id；理由只能依据该片段的摘要和问题原文，不要把摘要没有写的内容
   归到片段名下，也不要声称问题里写了它没写的东西。

严格按以下格式输出，每行一个片段，冒号后写一句适用理由：
<selected_knowledge>
片段id: 理由
</selected_knowledge>
"""
        else:
            return f"""Before an analysis is planned, you choose the fragments of the project's knowledge base that
apply to the user's question.
{context_block}
Candidate fragments (id | type | severity | lexical prefilter score | matched terms; the summary follows):
{catalogue}

Severity: must = ignoring it produces a wrong number; should = ignoring it produces a misleading
comparison or interpretation; info = background for interpreting the result.

Selection rules:
1. Choose the fragments that may change how this question is computed or how its result must be
   stated: the fragment constrains a computation step (units, weighting, a definition, the grid, the
   order of time averaging), or the answer has to state the caveat it raises.
2. Judge applicability from the summary, not from word overlap. The prefilter score is a hint only:
   a question mentioning a word does not make the fragment about that word apply.
3. The aim of this step is not to miss: the next step reads the full text of what you choose and
   drops what does not apply. When unsure whether a must or should fragment applies, include it.
   Rank generic points that every data read runs into below fragments about this question's own
   variables, datasets and methods.
4. Choose at most {max_selected}. If none applies, output the tags empty.
5. Invent nothing: choose only ids from the candidate list above, and ground each reason in that
   fragment's summary and the question's own words. Do not attribute to a fragment anything its
   summary does not say, and do not claim the question says something it does not.

Follow this format strictly, one fragment per line, with a one-sentence reason after the colon:
<selected_knowledge>
fragment-id: reason
</selected_knowledge>
"""

    @staticmethod
    def get_knowledge_selection_input(language: str, user_query: str) -> str:
        if language == 'zh':
            return f"用户问题：{user_query}\n\n现在输出选中的知识片段。"
        else:
            return f"User question: {user_query}\n\nNow output the fragments you choose."

    @staticmethod
    def get_knowledge_condense_prompt(language: str, sources: str, context_summary: str,
                                      max_notes: int, max_chars: int) -> str:
        """
        Build the system prompt that condenses the selected fragments into notes.

        sources is the full text of each selected fragment, headed by its id. The
        planner sees only the notes, so they carry what to do, not why. Numbers and
        `names` are to be quoted: src/docs/knowledge_retrieval.py checks every one
        against the fragment and the question, and a fragment whose notes fail falls
        back to its summary.
        """
        context_block = f"\n{context_summary}\n" if context_summary else ""

        if language == 'zh':
            return f"""你负责把知识片段压缩成简短的注意事项，交给数据分析的规划步骤。规划步骤只看得到你写的注意事项，看不到片段原文。
{context_block}
为本问题选出的片段（每段以 [片段id] 开头）：
<fragments>
{sources}
</fragments>

规则：
1. 每个片段只保留会改变【这道题】的计算方式或结论表述的内容。每个片段最多 {max_notes} 条，每条一行，
   最重要的一条写在最前面。
2. 每条都要能直接执行：写出要用的公式、阈值、单位换算、加权方式、工具参数或检查方法，或者答案必须
   声明的注意事项。不写背景、不写原文的论证，也不写本题不涉及的变量、模式或步骤。
3. 不准编造，内容必须服从片段原文：每条注意事项都要能在片段原文中找到依据；不要添加片段没有说的数字、
   工具参数、规则或例外，也不要写与原文相反的结论。数字、公式、变量名和工具名按片段原文照抄，变量名和
   工具名放在反引号里。可以引用问题本身定义的量（问题里的阈值、时段、区域）。
4. 只有当问题原文确实写出了与片段不同的定义或步骤时，才写一条说明：按问题的做法执行，并在答案中说明
   差异；这条说明要引用问题的原话。问题没有写到的地方，按片段执行，不要声称"问题规定了"什么。
5. 读完全文后发现某个片段不适用于本题，写 "[片段id] 不适用：原因"，不要为它写注意事项。
6. 全部输出控制在约 {max_chars} 个字符以内。

严格按以下格式输出：
<knowledge_notes>
[片段id]
- 注意事项
- 注意事项
[片段id] 不适用：原因
</knowledge_notes>
"""
        else:
            return f"""You condense knowledge fragments into short notes for the planning step of a data analysis. The
planner sees only your notes, not the fragments.
{context_block}
Fragments selected for this question (each headed by its [fragment-id]):
<fragments>
{sources}
</fragments>

Rules:
1. For each fragment keep only what changes how THIS question is computed or how its result must be
   stated. Write at most {max_notes} notes per fragment, one line each, the most important first.
2. Each note is actionable: the formula, threshold, unit conversion, weighting, tool argument or check
   to apply, or the caveat the answer must state. Leave out background, the source's reasoning, and
   points about variables, models or steps this question does not involve.
3. Invent nothing; the notes must follow the fragment text. Every note must be backed by what the
   fragment says: add no number, tool argument, rule or exception the fragment does not contain, and
   never state the opposite of what it says. Copy numbers, formulas, variable names and tool names
   exactly as the fragment writes them, and put variable and tool names in backticks. You may refer
   to quantities the question itself defines (its thresholds, periods, regions).
4. Only where the question's text actually states a definition or procedure that differs from the
   fragment, write a note saying so, quoting the question's words: the question's version is
   followed and the difference is stated in the answer. Where the question is silent, follow the
   fragment, and never claim the question prescribes something it does not say.
5. If, having read it, a fragment does not apply to this question, write
   "[fragment-id] not applicable: <reason>" and no notes for it.
6. Keep the whole answer under about {max_chars} characters.

Follow this format strictly:
<knowledge_notes>
[fragment-id]
- note
- note
[fragment-id] not applicable: reason
</knowledge_notes>
"""

    @staticmethod
    def get_knowledge_condense_input(language: str, user_query: str) -> str:
        if language == 'zh':
            return f"用户问题：{user_query}\n\n现在输出压缩后的注意事项。"
        else:
            return f"User question: {user_query}\n\nNow write the notes."

    @staticmethod
    def get_code_generation_prompt(language: str, full_specs: Dict, user_query: str,
                                   dataset_specs: Dict, conventions: Dict, plan: str = "") -> str:
        """
        Build the code-generation system prompt.
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
- 让错误抛出来。不要用 try/except 把某一步的异常记下来然后继续：那样进程会以「成功」结束，而答案里
  静静地缺了一块，修复环节也永远看不到这个问题。某个模式或某个变量按你写的方式处理不了时，正确做法
  是改调用，而不是捕获异常。只有当题目本身要求容忍缺失输入时才捕获，并且仍要把失败留在一个变量里，
  让最终答案能如实报告。
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
- 多子图共用色标时，色标上下限取汇总数据的稳健分位数，例如
  `vmax = float(np.nanpercentile(np.abs(np.concatenate([...])), 98))`，绝不要取绝对值最大值。
  近岸或边缘海常有比信号大一个量级的单点，`vmax = nanmax(abs(...))` 会把整个色标让给那一个格点，
  于是每个子图都渲染成一片空白。差值/距平用以零为中心的发散色图，量值用顺序色图；若做了截断，
  用 `extend='both'`（或 'max'/'min'）让色标两端显示箭头，并在图标题里说明截断方式。
- 色标标签只写物理量和单位，例如 "SSS bias (psu)"，控制在 40 个字符以内。竖直色标的标签比色标本身长时，
  保存的图片里两端会被裁掉。截断分位数、"共用色标"、数据来源等说明写进总标题或子图标题，不要写进色标标签。
- 给多个子图共用的 colorbar (`fig.colorbar(..., ax=axes)`) 之后，不要再调用 tight_layout()：
  它不知道这个 colorbar 的存在，会把子图铺到 colorbar 上。要么用
  `plt.subplots(..., constrained_layout=True)` 且完全不调 tight_layout，要么在添加 colorbar 之前调。
- 标注不要压住数据：图例用 `loc="best"` 或放到坐标区外。把图例钉在曲线峰值所在的角落，
  正好遮住这个子图要展示的特征。
- 曲线网格上的场，其二维坐标并不叫 lat/lon。请用 `da.attrs.get("lat_coord", "lat")` /
  `da.attrs.get("lon_coord", "lon")` 读出坐标名，不要直接索引 `da.coords["lat"]`——
  在 CESM2 和 MPI-ESM1-2-LR 上会抛 KeyError。
- 画曲线网格上的场（CESM2、MPI-ESM1-2-LR 的原生输出）的地图时，先用
  `regrid_horizontal(field, target="woa23")` 放到规则网格上（仅用于显示，所有汇报的数值仍在原生网格上算），
  再用一维 lon/lat 画 `pcolormesh`。不要把它的二维坐标直接交给 `pcolormesh`/`contourf`：模式经度在网格
  接缝处回绕跳变，会把场抹成横向条带并画错位置。
- 画跨日界线的区域（如北太平洋 110E-100W、赤道太平洋 120E-280E）时，用连续的 0-360 经度轴：仅为作图，
  对规则网格场的副本做 `p = da.assign_coords(lon=da.lon % 360).sortby("lon")`，再把 x 轴范围设成该区域。
  在 -180..180 的经度轴上，这个区域会被劈成子图两端的两条，中间一大片空白。
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
- Let failures raise. Do NOT wrap a step in try/except to record the error and carry on: the run
  then exits cleanly with part of the answer silently missing, and the repair pass never sees the
  problem. If one model or one variable cannot be processed the way you wrote it, the right
  response is to fix the call, not to catch the exception. Only catch one when the question itself
  asks you to tolerate a missing input, and then still leave the failure in a variable the answer
  can report.
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
- Colour limits on a shared scale: take them from a robust percentile of the pooled data, e.g.
  `vmax = float(np.nanpercentile(np.abs(np.concatenate([...])), 98))`, and never from the maximum
  absolute value. A single coastal or marginal-sea cell is routinely an order of magnitude larger
  than the signal, and `vmax = nanmax(abs(...))` then hands the whole scale to that one cell so
  every panel renders as uniform white. Use a diverging map centred on zero for a difference or
  anomaly and a sequential map for a magnitude. When values are clipped, pass `extend='both'` (or
  'max'/'min') so the colourbar shows arrows, and say how they were clipped in the figure title.
- Keep a colourbar label to the quantity and its units, e.g. "SSS bias (psu)", under about 40
  characters. A vertical colourbar label longer than the bar is cut off at both ends in the saved
  image. Notes such as the clipping percentile, "shared scale" or the data source go in the figure
  or panel title, never in the colourbar label.
- Do NOT call tight_layout() after attaching a colorbar to several axes (`fig.colorbar(..., ax=axes)`)
  — it does not know about that colorbar and will lay the panels over it. Either build the figure
  with `plt.subplots(..., constrained_layout=True)` and add no tight_layout at all, or call
  tight_layout BEFORE adding the colorbar.
- Keep annotation off the data: `loc="best"` for a legend, or place it outside the axes. A legend
  pinned to a corner where a curve peaks hides the very feature the panel is there to show.
- A field on a curvilinear grid has 2-D coordinates that are NOT called lat/lon. Read the names off
  the array with `da.attrs.get("lat_coord", "lat")` / `da.attrs.get("lon_coord", "lon")` rather than
  indexing `da.coords["lat"]`, which raises KeyError on CESM2 and MPI-ESM1-2-LR.
- To map a field on a curvilinear grid (CESM2 or MPI-ESM1-2-LR native output), first put it on a
  regular grid with `regrid_horizontal(field, target="woa23")` (for display only; every reported
  number stays on the native grid) and draw it with `pcolormesh` on the 1-D lon/lat. Do not hand its
  2-D coordinates to `pcolormesh`/`contourf`: the model longitude wraps and jumps at the grid seam,
  which smears the field into horizontal bands and puts it in the wrong place.
- To map a region that crosses the dateline (the North Pacific 110E-100W, the equatorial Pacific
  120E-280E), use a continuous 0-360 longitude axis: for the plot only, on a copy of the regular-grid
  field, `p = da.assign_coords(lon=da.lon % 360).sortby("lon")`, then set the x limits to the region.
  On a -180..180 axis the region is split into two strips at the panel's edges with empty space
  between them.
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
    def get_fix_prompt(language: str, code: str, exec_result: Dict, user_query: str,
                       full_specs: Dict, code_type: str,
                       upstream_code: str = "", available_variables: Dict = None) -> str:
        """
        Build the repair prompt for one failed code block.
        """
        upstream_block = ""
        variables_block = ""
        if upstream_code:
            upstream_block = (
                f"""
已经成功执行的上游代码（它定义的变量在当前命名空间中已经可用）：
<upstream_code>
{upstream_code}
</upstream_code>
""" if language == 'zh' else f"""
The code that already ran successfully before this block (the variables it defines are
already in the namespace):
<upstream_code>
{upstream_code}
</upstream_code>
""")
        if available_variables:
            variables_block = (
                f"""
当前命名空间中已有的变量（名称 -> 内容摘要），可以直接使用：
{json.dumps(available_variables, indent=2, ensure_ascii=False, default=str)}
""" if language == 'zh' else f"""
Variables already present in the namespace (name -> summary of its contents); use them
as they are:
{json.dumps(available_variables, indent=2, ensure_ascii=False, default=str)}
""")

        # Only rules the prompt can back up: a reuse rule with no list of variables
        # under it reads as an instruction to invent one.
        if language == 'zh':
            focus = ["修复上述具体错误"]
            if available_variables:
                focus.append("变量名必须取自上面的可用变量清单；清单里已有的量直接引用，"
                             "不要重新加载数据或重新计算")
            else:
                focus.append("确保变量名与函数调用正确")
            if upstream_code:
                focus.append("不要输出上游代码，也不要修改它")
            focus.append("使用正确的语法、函数签名和库函数")
        else:
            focus = ["Fixing the specific error mentioned above"]
            if available_variables:
                focus.append("Taking variable names from the list of available variables above; "
                             "use what is already there instead of reloading the data or "
                             "recomputing a quantity that exists")
            else:
                focus.append("Ensuring proper variable names and function calls")
            if upstream_code:
                focus.append("Not reproducing or modifying the upstream code")
            focus.append("Using correct syntax, function signatures and library functions")
        focus_block = "\n".join(f"{i}. {line}" for i, line in enumerate(focus, 1))

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
{upstream_block}{variables_block}
可用工具及其规格说明：
{json.dumps(full_specs, indent=2, ensure_ascii=False)}

请只修复上面这一个 {code_type} 代码块，并在 <execute></execute> 标签内输出修正后的完整版本。
重点关注：
{focus_block}
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
{upstream_block}{variables_block}
Available tools and their specifications:
{json.dumps(full_specs, indent=2)}

Fix this one {code_type} block only, and output the corrected, complete version inside
<execute></execute> tags.
Focus on:
{focus_block}
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
