# JingweiAnalyst

JingweiAnalyst answers ocean and climate analysis questions by generating and executing Python code against local datasets. The full system combines expertise tools with a curated knowledge base, then saves the answer, figures, code, and execution records. This repository also includes three comparison modes, ten benchmark tasks, reference answers, and evaluation scripts.

To get started, prepare the data, install the environment, configure a model, and run a small task batch. Run all commands below from the repository root.

## 1. Prepare the data

Download the datasets needed by tasks and arrange them as below. Paths in the table are relative to your chosen data storage directory, following the archive's `ocean_data/` and `benchmark_cmip6/` layout. The system reads local files and does not download them automatically.

| Dataset ID | Download source | Required local organization |
| --- | --- | --- |
| `woa23` | NOAA [WOA23 data](https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/) and [masks](https://www.ncei.noaa.gov/data/oceans/woa/WOA23/MASKS/) | Set `root` to `ocean_data/WOA23/climatology/`. Store each file under `<variable>/<period-subdirectory>/woa23_<tag>_<prefix><nn>_01.nc`, using the exact subdirectories and tags in `datasets.yaml`. Set `masks_root` to `ocean_data/WOA23/masks/`, containing `basinmask_01.msk` and `landsea_01.msk`. |
| `glodap` | [GLODAPv2.2016b mapped archive](https://glodap.info/index.php/mapped-data-product/) | Extract the mapped-climatology archive and set `root` to `ocean_data/GLODAPv2.2016b_MappedClimatologies/`. Place `GLODAPv2.2016b.<varfile>.nc` directly in this directory, e.g. `GLODAPv2.2016b.TCO2.nc`. Do not add per-variable subdirectories. |
| `seanoe_mld` | SEANOE [v2024](https://doi.org/10.17882/98226) and [v2023](https://doi.org/10.17882/91774) | Set `root` to `ocean_data/SEANOE/`. Place `mld_dr003_ref5m_v2024.nc` and `mld_dr003_ref10m_v2023.nc` directly inside it. The configured default is v2024; v2023 is needed when explicitly requested. |
| `ersst` | NOAA PSL [sst.mnmean.nc](https://psl.noaa.gov/thredds/fileServer/Datasets/noaa.ersst.v5/sst.mnmean.nc) | Set `root` to `ocean_data/ERSSTv5/`, containing the single combined monthly file `sst.mnmean.nc`. Individual monthly files are not the input layout used by this reader. |
| `cpc_indices` | NOAA CPC [1991-2020 Nino](https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii), [1981-2010 Nino](https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.81-10.ascii), [ONI](https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt), [OI SST](https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices); NOAA NCEI [PDO directory](https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/v6/index/) | Set `root` to `ocean_data/CPC/`. Keep `ersst5.nino.mth.91-20.ascii`, `ersst5.nino.mth.81-10.ascii`, `oni.ascii.txt`, `sstoi.indices`, and `ersst.v6.pdo.dat` directly inside it, preserving filenames and text formats. |
| `cmip6_historical` | [CMIP / ESGF data access](https://www.wcrp-cmip.org/cmip-data-access/) | Set `root` to `benchmark_cmip6/historical/`. Place NetCDF files directly under `<model>/<variable>/`, including fixed fields in their own variable directories. |
| `cmip6_ssp245` | [CMIP / ESGF data access](https://www.wcrp-cmip.org/cmip-data-access/) | Set `root` to `benchmark_cmip6/ssp245/`, with the same `<model>/<variable>/*.nc` layout. Retain historical fixed fields for fallback where the scenario's variable directory is absent. |

### WOA23 and CMIP6 directory details

WOA23 uses the configured period subdirectories, not the directory names on the download server. For example:

```text
ocean_data/WOA23/
  climatology/
    temperature/1971-2000_climate_normal/woa23_decav71A0_t00_01.nc
    temperature/1991-2020_climate_normal/woa23_decav91C0_t01_01.nc
    salinity/1971-2000_climate_normal/woa23_decav71A0_s00_01.nc
    oxygen/all_time_periods/woa23_all_o00_01.nc
    nitrate/all_time_periods/woa23_all_n00_01.nc
  masks/
    basinmask_01.msk
    landsea_01.msk
```

Use `nn=00` for annual, `01`-`12` for monthly, and `13`-`16` for seasonal files. The registered 1991-2020 temperature and salinity directories contain monthly files only. Nutrients use `all_time_periods`; the registry lists every supported variable/period combination. Preserve case in filenames, including the `A` prefix for apparent oxygen utilization and `O` for oxygen saturation. Keep `nitrate/all_time_periods/woa23_all_n00_01.nc` available when using tools that read the shared WOA depth reference.

For CMIP6, use model directory names `CESM2`, `GFDL-ESM4`, and `MPI-ESM1-2-LR`, with grid labels `gn`, `gr`, and `gn`, respectively. For example:

```text
benchmark_cmip6/
  historical/
    CESM2/
      tos/tos_Omon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc
      areacello/areacello_Ofx_CESM2_historical_r1i1p1f1_gn.nc
      sftof/sftof_Ofx_CESM2_historical_r1i1p1f1_gn.nc
    GFDL-ESM4/
      tos/tos_Omon_GFDL-ESM4_historical_r1i1p1f1_gr_185001-186912.nc
      tos/...
    MPI-ESM1-2-LR/
      tos/tos_Omon_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001-186912.nc
      tos/...
  ssp245/
    CESM2/
      thetao/...
    GFDL-ESM4/
      thetao/...
    MPI-ESM1-2-LR/
      thetao/...
```

The examples show the directory structure, not a complete download list. Put all time chunks needed for a task in the corresponding variable directory. The reader searches only `<root>/<model>/<variable>/*.nc`: remove intermediate ESGF table/member/grid/version directories when organizing downloads. Keep original filenames, one ensemble member, one grid, and one data version per variable directory; mixed ensemble members are rejected.

Include `areacello` for model area weighting and other grid fields needed by the task (`deptho`, `sftof`, and `volcello`, where available). CESM2 also needs `sftof` when its fixed fields are masked. GFDL-ESM4 `volcello` is a time-dependent `Omon` field, so retain its time chunks. For SSP2-4.5, fallback to historical fixed fields applies only when the corresponding scenario variable directory does not exist; an empty directory does not trigger fallback.

### Set paths in the configuration

In [configs/datasets.yaml](configs/datasets.yaml), replace every dataset's `root` and the separate `woa23.masks_root` with the corresponding absolute paths on your machine. For example, edit these fields in the existing registry rather than replacing the full file:

```yaml
datasets:
  woa23:
    root: /path/to/data/ocean_data/WOA23/climatology
    masks_root: /path/to/data/ocean_data/WOA23/masks
  glodap:
    root: /path/to/data/ocean_data/GLODAPv2.2016b_MappedClimatologies
  seanoe_mld:
    root: /path/to/data/ocean_data/SEANOE
  ersst:
    root: /path/to/data/ocean_data/ERSSTv5
  cpc_indices:
    root: /path/to/data/ocean_data/CPC
  cmip6_historical:
    root: /path/to/data/benchmark_cmip6/historical
  cmip6_ssp245:
    root: /path/to/data/benchmark_cmip6/ssp245
```

Preserve the registry's filenames, subdirectories, and metadata when using this layout. Its CMIP6 `layout` and `filename_template` fields describe the archive; changing them does not change the reader's directory search. Consult [cases/tasks](cases/tasks) for each task's variables and periods, and the registry notes for missing variables and ensemble-member exceptions. Different data versions or ensemble members can change results relative to the reference answers.

The release archive includes `assets/` and `cache/`. After extracting it, point [configs/default_paths.yaml](configs/default_paths.yaml) to those directories and ensure the cache is writable:

```yaml
cache_dir: /absolute/path/to/JingweiAnalyst/cache
natural_earth_path: /absolute/path/to/JingweiAnalyst/assets/NaturalEarth
```

## 2. Set up the environment

The dependency baseline is Python 3.11.14 on Linux. Install Conda first, then create the environment:

```bash
conda create -n jingweianalyst -c conda-forge python=3.11.14 pip esmpy=8.6.1
conda activate jingweianalyst
python -m pip install -r requirement.txt
python -m pip check
python -c "import esmpy, xesmf, netCDF4, geopandas"
```

[requirement.txt](requirement.txt) pins the direct dependencies used by the project. It is not a complete environment lock file. ESMPy supplies the native regridding backend for xESMF and must be installed separately through Conda or from source; see the [xESMF installation guide](https://xesmf.readthedocs.io/en/stable/installation.html).

The default model interface supports OpenAI-compatible APIs, including DeepSeek, OpenAI, and Qwen. For other providers, install the corresponding integration:

```bash
# Only needed for Anthropic / Claude
python -m pip install -r requirement.txt langchain-anthropic

# Only needed for Google / Gemini
python -m pip install -r requirement.txt langchain-google-genai
```

## 3. Configure and run the system

### Shared settings and a first task

Edit [configs/agent_full.yaml](configs/agent_full.yaml). The other three mode configurations inherit its settings through `extends`.

1. Set `llm.provider` and `llm.model_name` to a model available through your account. Set `base_url` if using a proxy or custom endpoint, and use a temperature supported by the model (`null` omits it).
2. Export the API key in the same shell used to run the system. For example:

   ```bash
   export DEEPSEEK_API_KEY="your-api-key"
   ```

   Other defaults are `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `GOOGLE_API_KEY`, and `DASHSCOPE_API_KEY`. Set `llm.api_key_env` to use a different variable name. Keep keys out of configuration files.
3. Set `agent.result_dir` in **each mode's configuration** to a new output directory. The supplied values point to archived experiment directories.
4. Set `batch_run.test_folder` and `batch_run.repeats`. The supplied defaults run all ten benchmark tasks in `cases/tasks` three times; use `repeats: 1` for a first run.

For a first run, use the existing [T1: multi-model ENSO amplitude assessment](cases/tasks/T1_cmip6_enso_amplitude_multimodel.json). It compares CESM2, GFDL-ESM4, and MPI-ESM1-2-LR historical `tos` against ERSSTv5 over January 1950 to December 2014. Prepare `tos` covering that period and `areacello` for all three models, CESM2 `sftof` for the area mask, and the ERSSTv5 file.

After configuring the full system and its output directory, run this example directly from the existing task file:

```bash
python - <<'PYTHON'
import json
from pathlib import Path
from src.models.agent_with_tools import AgentWithTools, TOOLS
from src.tools.data_openers import list_datasets

path = Path("cases/tasks/T1_cmip6_enso_amplitude_multimodel.json")
task = json.loads(path.read_text(encoding="utf-8"))[0]
agent = AgentWithTools(TOOLS, list_datasets(),
                       config_path=Path("configs/agent_full.yaml"))
session_id = agent.create_session(label=task["case_id"])
try:
    agent.run(task["query"], session_id=session_id)
finally:
    agent.end_session(session_id)
PYTHON
```

This runs T1 once using the full configuration, independently of `batch_run.repeats`. The commands below use the batch runner: `batch_run.test_folder` must be a **directory**, not a JSON file path. With `test_folder: cases/tasks`, they run all supplied tasks. Each JSON file contains a list of questions; each item starts a new session.

### Full: tools and documents

The full system retrieves knowledge, plans the analysis, selects scientific tools, generates and executes code, repairs failed code within its retry budget, and writes an answer.

In [configs/agent_full.yaml](configs/agent_full.yaml), keep `arm: full` and `knowledge.enabled: true`; set `agent.result_dir`, for example, to `results/my_full`.

```bash
python -m src.models.agent_with_tools --config configs/agent_full.yaml
```

### +Tools: tools without documents

This mode retains the expertise tools but disables knowledge retrieval. It corresponds to the `no_knowledge` arm (also called `-K`).

In [configs/agent_with_tools_no_knowledge.yaml](configs/agent_with_tools_no_knowledge.yaml), keep `arm: no_knowledge` and `knowledge.enabled: false`; set `agent.result_dir`, for example, to `results/my_tools`.

```bash
python -m src.models.agent_with_tools --config configs/agent_with_tools_no_knowledge.yaml
```

### +Docs: documents without project tools

This mode retains knowledge retrieval. Generated code reads and processes the raw data with installed Python libraries instead of the project's scientific tools. It corresponds to the `no_tools` arm (also called `-T`).

In [configs/agent_with_knowledge_no_tools.yaml](configs/agent_with_knowledge_no_tools.yaml), keep `arm: no_tools`; knowledge retrieval is inherited from the full configuration. Set `agent.result_dir`, for example, to `results/my_docs`.

```bash
python -m src.models.ablation_tools.agent_without_tools --config configs/agent_with_knowledge_no_tools.yaml
```

### Direct LLM

This baseline generates code directly from the question and dataset descriptions, without a separate planning step, knowledge retrieval, or project tools. It still executes code, attempts repairs, and produces an answer.

In [configs/direct_llm.yaml](configs/direct_llm.yaml), keep `arm: direct_llm` and `knowledge.enabled: false`; set `agent.result_dir`, for example, to `results/my_direct`.

```bash
python -m src.models.direct_llm.agent_direct_llm --config configs/direct_llm.yaml
```

| Mode | Configuration `arm` | Expertise tools | Knowledge retrieval |
| --- | --- | --- | --- |
| Full | `full` | Yes | Yes |
| +Tools | `no_knowledge` | Yes | No |
| +Docs | `no_tools` | No | Yes |
| Direct LLM | `direct_llm` | No | No |

Configuration mappings merge recursively; other values replace inherited values. Shared changes belong in `agent_full.yaml`, while mode-specific overrides belong in the child files. Relative runtime paths resolve against the repository root; `extends` resolves relative to its configuration file.

For comparisons, Full versus +Tools measures the contribution of retrieval, and Full versus +Docs measures the contribution of project tools.

## 4. Inspect results and evaluate

Each session writes turn directories under `agent.result_dir`, such as `<case>_<model>_<timestamp>/turn01/`. Start with `solution.md` and `figures/`. The same directory records:

- `code/compute.py`, `code/visualize.py`, and `code/attempts/`: generated code and repair history.
- `observation.json`: execution results summarized for the model.
- `prompts/`, `plan.md`, and `knowledge.json`: model inputs, planning, and retrieval records where applicable.
- `agent_config.json` and `manifest.json`: configuration, provenance, and execution status.

Failed runs may contain only the files written before failure. Check the execution records as well as the final answer.

To rerun saved computation and plotting code with the same local data and environment:

```bash
python results/my_full/<session>/turn01/code/run.py
```

Replace `<session>` with an actual session directory. The system executes generated Python in the current environment; use an environment where that code is permitted to access your local data.

### Benchmark answer evaluation

Configure the extraction model in [configs/evaluate.yaml](configs/evaluate.yaml), then score a benchmark output directory:

```bash
python -m src.eval.evaluate --config configs/evaluate.yaml --result-dir results/my_full
```

The evaluator uses an LLM to extract reported values and compares them with [cases/reference_ans](cases/reference_ans) using `grading_spec.yaml`. It requires the extractor's API key and writes `evaluation.json` by default. Custom questions need corresponding reference answers and grading entries before they can be scored.

For evaluation figures, update the result directories and input filenames in [configs/visualize.yaml](configs/visualize.yaml), then run:

```bash
python -m src.eval.visualize --config configs/visualize.yaml
```

The supplied figure configuration targets archived experiments; some panels also require figure/report assessments or loss-analysis files in addition to numerical evaluation results.

### Scientific tool verification

The verification suite checks analytic solutions, independent implementations, published products, physical magnitudes, consistency, Monte Carlo properties, and safeguards. It uses local datasets and does not require an LLM API key.

```bash
# Run all categories sequentially in separate processes
python -u -m src.eval.tool_verification.run_all --out results/tool_verification/my_run

# Run one selected tool within a category
python -u -m src.eval.tool_verification.run_category cat2_reference --out results/tool_verification/my_subset --only compute_pdo_index

# Rebuild the report from category records
python -m src.eval.tool_verification.summarize results/tool_verification/my_run
```

Outputs include per-category records and logs, `report.md`, and `summary.json`. The sequential runner avoids concurrent NetCDF reads across categories, which caused native crashes on the original NFS storage. Inspect reported failures and errors rather than relying only on process completion.

### Knowledge base and repository layout

The curated documents live in [docs/knowledge](docs/knowledge). After editing them, validate the catalogue:

```bash
python -m src.docs.knowledge_index --check
```

Run the same command without `--check` to validate and regenerate `INDEX.md`. Full and +Docs validate the knowledge base at startup by default.

| Directory | Purpose |
| --- | --- |
| `configs/` | Dataset registry, resource paths, model settings, and evaluation settings |
| `src/models/` | Agents, prompts, and model-provider interfaces |
| `src/tools/`, `src/utils/` | Scientific operations, data readers, and geographic utilities |
| `src/docs/`, `docs/knowledge/` | Knowledge indexing, retrieval, and source documents |
| `cases/tasks/`, `cases/reference_ans/` | Benchmark questions, reference calculations, and grading specifications |
| `src/eval/` | Answer evaluation, plotting, and scientific tool verification |
| `results/` | Run artifacts and evaluation outputs |

When reporting results, retain the run configuration and data product versions, and cite the original data providers. Model responses and updated data downloads can differ from the archived experiments.
