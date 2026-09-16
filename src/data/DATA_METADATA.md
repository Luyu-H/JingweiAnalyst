# 数据集元数据总览（Data Metadata Reference）

> 生成日期：2026-09-09
> 覆盖根目录：`/data1/lorraine-20419/`
> 元数据提取方式：`xarray` / `netCDF4` 直接读取文件头（未做数值统计），ASCII 文件按文本解析。
> 配套文件：`cmip6_metadata_raw.json` —— §1/§2 中每个模式-情景-变量的机读元数据（维度、坐标范围与步长、变量属性、全局属性、文件清单与体积）。

## 0. 总目录一览

| # | 路径 | 类型 | 体积 | 文件格式 |
|---|------|------|------|----------|
| 1 | `benchmark_cmip6/historical` | CMIP6 模式历史模拟 | 517 GB | NetCDF4 |
| 2 | `benchmark_cmip6/ssp245` | CMIP6 SSP2-4.5 情景 | 291 GB | NetCDF4 |
| 3 | `ocean_data/CPC` | NOAA/CPC ENSO 指数 | 228 KB | ASCII |
| 4 | `ocean_data/ERSSTv5` | 重建月平均 SST | 152 MB | NetCDF |
| 5 | `ocean_data/GLODAPv2.2016b_MappedClimatologies` | 全球海洋碳化学气候态 | 1.4 GB | NetCDF |
| 6 | `ocean_data/SEANOE` | 混合层深度气候态 | 50 MB | NetCDF + PDF |
| 7 | `ocean_data/WOA23` | 世界大洋图集 2023 | 3.5 GB | NetCDF + ASCII mask |

**共同要点**
- 所有格点数据经纬度均为地理坐标（WGS84 / EPSG:4326 等价）；但**经度约定不统一**：CMIP6 与 ERSSTv5 为 0–360°E，WOA23 与 SEANOE 为 −180–180°E，GLODAP 为 20.5–379.5°E（以 20°E 为起点绕行一周）。跨数据集比对前必须先统一经度。
- 深度坐标一律 `positive: down`，但**单位不统一**：CESM2 为 **厘米**，其余（GFDL、MPI、WOA23、GLODAP）为 **米**。
- 时间轴日历不统一：CESM2 = `noleap`/`365_day`，GFDL-ESM4 = `noleap`，MPI-ESM1-2-LR = `proleptic_gregorian`。做时间对齐时须显式处理。

---

# 1. `/data1/lorraine-20419/benchmark_cmip6/historical`

## 1.1 数据集概况

CMIP6 `historical` 试验（全强迫历史模拟），海洋分量 `Omon`（月频）与 `Ofx`（固定场）。

| 模式 | 机构 | 网格标签 | 网格描述 | 名义分辨率 | 空间维数 | 垂向层数 | 体积 |
|------|------|----------|----------|-----------|---------|---------|------|
| CESM2 | NCAR | `gn` | 原生 gx1v7 位移极点网格 (384×320) | 100 km | nlat=384, nlon=320 | 60 (POP2) | 303 GB |
| GFDL-ESM4 | NOAA-GFDL | `gr` | 由 tripolar（名义 0.5°）保守插值到规则 1°×1° (180×360) | 1×1 degree | lat=180, lon=360 | 35 | 87 GB |
| MPI-ESM1-2-LR | MPI-M | `gn` | 原生曲线网格 (220×256) | 250 km | j=220, i=256 | 40 | 127 GB |

- Conventions：CESM2 / MPI 为 `CF-1.7 CMIP-6.2`；GFDL-ESM4 为 `CF-1.7 CMIP-6.0 UGRID-1.0`。
- 文件命名：`<var>_<table>_<model>_<experiment>_<variant>_<grid>[_<YYYYMM-YYYYMM>].nc`

## 1.2 时间维

| 模式 | 时间范围 | 时次数 | 日历 | 分块 |
|------|---------|-------|------|------|
| CESM2 | 1850-01-15 → 2014-12-15 | 1980 月 | `noleap` | 单文件（1 个 .nc） |
| GFDL-ESM4 | 1850-01-16 → 2014-12-16 | 1980 月 | `noleap` | 9 个文件，每 20 年（末段 2010-2014 为 5 年） |
| MPI-ESM1-2-LR | 1850-01-16 → 2014-12-16 | 1980 月 | `proleptic_gregorian` | 9 个文件，每 20 年（末段 5 年） |

时间戳为**月中心时刻**，附 `time_bnds`（月起止边界）。

## 1.3 空间与垂向坐标

**CESM2**（曲线网格，需用二维 `lat`/`lon` 定位）
- `nlat`(384)、`nlon`(320) 为**索引**（1..384 / 1..320），非地理坐标；
- 二维 `lat(nlat,nlon)` ∈ [−79.22, 89.71] °N，`lon(nlat,nlon)` ∈ [0.015, 360.0] °E；带 4 顶点 `lat_bnds`/`lon_bnds`；
- `lev`(60)：**单位 centimeters**，500 → 537500 cm（即 5 m → 5375 m），上层 10 m 等间距（Δ=1000 cm），深层最大 Δ≈24997 cm。
  首层：5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155 m …
- 例外：`chl` 仅 15 层（500–14500 cm，即 5–145 m，仅上 150 m）。

**GFDL-ESM4**（规则经纬网格，最易用）
- `lat`(180) = −89.5 … 89.5，Δ=1.0°（均匀）；
- `lon`(360) = 0.5 … 359.5，Δ=1.0°（均匀）；
- `lev`(35)：单位 **m**，2.5 → 6500 m，非均匀（Δ 7.5 → 500 m）：
  `2.5, 10, 20, 32.5, 51.25, 75, 100, 125, 156.25, 200, 250, 312.5, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1537.5, 1750, 2062.5, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500`

**MPI-ESM1-2-LR**（曲线网格）
- `j`(220)、`i`(256) 为索引（0..219 / 0..255）；
- 二维 `latitude(j,i)` ∈ [−83.97, 89.73] °N，`longitude(j,i)` ∈ [0.007, 360.0] °E；带 4 顶点边界；
- `lev`(40)：单位 **m**，6 → 5720 m，非均匀（Δ 10 → 550 m）：
  `6, 17, 27, 37, 47, 57, 68.5, 82.5, 100, 122.5, 150, 182.5, 220, 262.5, 310, 362.5, 420, 485, 560, 645, 740, 845, 960, 1085, 1220, 1365, 1525, 1700, 1885, 2080, 2290, 2525, 2785, 3070, 3395, 3770, 4195, 4670, 5170, 5720`
- **注意**：`wo`、`wmo` 使用 41 层（`lev_2`，界面层），非 40。

## 1.4 变量清单（historical）

`✓` = 存在；`—` = 该模式目录下无此变量。

| 变量 | 长名 | 单位 | 维度 | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|------|------|------|------|:---:|:---:|:---:|
| `thetao` | Sea Water Potential Temperature | degC | (t,z,y,x) | ✓ | ✓ | ✓ |
| `so` | Sea Water Salinity | 0.001（≈ psu） | (t,z,y,x) | ✓ | ✓ | ✓ |
| `tos` | Sea Surface Temperature | degC | (t,y,x) | ✓ | ✓ | ✓ |
| `mlotst` | Ocean Mixed Layer Thickness Defined by Sigma T | m | (t,y,x) | ✓ | ✓ | ✓ |
| `uo` | Sea Water X Velocity | m s-1 | (t,z,y,x) | ✓ | — | ✓ |
| `vo` | Sea Water Y Velocity | m s-1 | (t,z,y,x) | ✓ | — | ✓ |
| `wo` | Sea Water Vertical Velocity | m s-1 | (t,z,y,x) | ✓ | — | ✓ (41 层) |
| `umo` | Ocean Mass X Transport | kg s-1 | (t,z,y,x) | — | — | ✓ |
| `vmo` | Ocean Mass Y Transport | kg s-1 | (t,z,y,x) | — | — | ✓ |
| `wmo` | Upward Ocean Mass Transport | kg s-1 | (t,z,y,x) | — | — | ✓ (41 层) |
| `dissic` | Dissolved Inorganic Carbon Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `dissoc` | Dissolved Organic Carbon Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `talk` | Total Alkalinity | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `ph` | pH | 1（无量纲） | (t,z,y,x) | ✓ | ✓ | ✓ |
| `o2` | Dissolved Oxygen Concentration | mol m-3 | (t,z,y,x) | — | ✓ | ✓ |
| `o2sat` | Dissolved Oxygen Concentration at Saturation | mol m-3 | (t,z,y,x) | — | ✓ | ✓ |
| `no3` | Dissolved Nitrate Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `po4` | Total Dissolved Inorganic Phosphorus Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `si` | Total Dissolved Inorganic Silicon Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `dfe` | Dissolved Iron Concentration | mol m-3 | (t,z,y,x) | ✓ | ✓ | ✓ |
| `chl` | Mass Concentration of Total Phytoplankton expressed as Chlorophyll | kg m-3 | (t,z,y,x) | ✓ (15 层) | ✓ | ✓ |
| `phycos` | (Sea) Surface Phytoplankton Carbon Concentration | mol m-3 | (t,y,x) | ✓ | ✓ | ✓ |
| `intdic` | Dissolved Inorganic Carbon Content（水柱积分） | kg m-2 | (t,y,x) | — | ✓ | ✓ |
| `intdoc` | Dissolved Organic Carbon Content | kg m-2 | (t,y,x) | ✓ | ✓ | ✓ |
| `intpoc` | Particulate Organic Carbon Content | kg m-2 | (t,y,x) | — | ✓ | ✓ |
| `areacello` | Grid-Cell Area for Ocean Variables | m2 | (y,x) | ✓ | ✓ | ✓ |
| `volcello` | Ocean Grid-Cell Volume | m3 | (z,y,x) | ✓ Ofx | ✓ **Omon**（含时间轴） | ✓ Ofx |
| `deptho` | Sea Floor Depth Below Geoid | m | (y,x) | ✓ | ✓ | ✓ |
| `sftof` | Sea Area Percentage | % | (y,x) | ✓ | — | ✓ |

## 1.5 需注意的不一致（historical）

1. **成员号（variant_label）不统一**：
   - `historical/CESM2/phycos` 为 **r6i1p1f1**，其余变量为 r1i1p1f1；
   - `historical/MPI-ESM1-2-LR/ph` 为 **r10i1p1f1**，其余变量为 r1i1p1f1。
   → 与其他变量做逐格点/逐时刻联合分析时，这两个变量并非同一集合成员，不可视为同一条模拟轨迹。
2. **`volcello` 在 GFDL-ESM4 中是 `Omon`（随时间变化，9 个分块文件，3.6 GB）**，而在 CESM2 / MPI 中是 `Ofx` 固定场。做体积加权积分时代码需分支处理。
3. CESM2 缺 `o2` / `o2sat` / `intdic` / `intpoc`；GFDL-ESM4 缺 `uo/vo/wo`、`umo/vmo/wmo`、`sftof`。
4. `so` 单位标注为 `0.001`（CMIP6 对实用盐度的标准写法），不是 psu 字符串。

---

# 2. `/data1/lorraine-20419/benchmark_cmip6/ssp245`

## 2.1 数据集概况

CMIP6 ScenarioMIP `ssp245`（SSP2-4.5，2100 年辐射强迫约 4.5 W/m²，浓度驱动）。网格、坐标定义与 §1 的 historical 完全一致（同一模式共用网格），差别在时间范围、成员号与变量覆盖。

| 模式 | 时间范围 | 时次数 | 日历 | 分块 | 体积 |
|------|---------|-------|------|------|------|
| CESM2 | 2015-01-15 → 2100-12-15 | 1032 月 | `365_day` | 2 个文件（2015-2064, 2065-2100） | 180 GB |
| GFDL-ESM4 | 2015-01-16 → 2100-12-16 | 1032 月 | `noleap` | 5 个文件，每 20 年（末段 2095-2100 为 6 年） | 45 GB |
| MPI-ESM1-2-LR | 2015-01-16 → 2100-12-16 | 1032 月 | `proleptic_gregorian` | 5 个文件，每 20 年（末段 6 年） | 66 GB |

## 2.2 变量清单（ssp245）

变量的长名 / 单位 / 维度与 §1.4 表完全相同，下表只标注覆盖情况差异。

| 变量 | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|------|:---:|:---:|:---:|
| `thetao` `so` `dissic` `dissoc` `talk` `ph` `o2` `o2sat` `no3` `po4` `si` `dfe` `chl` `phycos` `intdic` `intdoc` | ✓ | ✓ | ✓ |
| `intpoc` | — | ✓ | ✓ |
| `uo` `vo` `wo` | ✓ | — | ✓ |
| `umo` `vmo` `wmo` | — | — | ✓ |
| `tos` | — | — | — |
| `mlotst` | — | — | — |
| `areacello` `deptho` `sftof` | ✓ | — | ✓ |
| `volcello` | ✓ (Ofx) | ✓ (**Omon**) | ✓ (Ofx) |

## 2.3 需注意的不一致（ssp245）

1. **`tos` 与 `mlotst` 在三个模式的 ssp245 目录下均缺失**（historical 有）。若需要情景期的 SST / 混合层深度，须另行下载，或用 `thetao` 最上层近似 `tos`。
2. **CESM2 的 ssp245 主体使用成员 `r10i1p1f1`**，而 `deptho`、`sftof` 两个固定场是 **r1i1p1f1**。固定场与网格无关成员差异，通常可用，但需知悉。
3. CESM2 的 ssp245 日历属性写作 `365_day`，historical 写作 `noleap`（二者等价，但字符串不同，硬编码判断时会踩坑）。
4. GFDL-ESM4 的 ssp245 目录**不含任何 `Ofx` 固定场**（`areacello`/`deptho`/`sftof` 均无），需借用 historical 目录下的同名文件（网格相同）。

---

# 3. `/data1/lorraine-20419/ocean_data/CPC`

NOAA/NWS 气候预测中心（CPC）ENSO 监测指数，纯 ASCII 表格，共 228 KB。

| 文件 | 内容 | 时间范围 | 行数 | 列 |
|------|------|---------|-----|-----|
| `oni.ascii.txt` | ONI（Oceanic Niño Index），Niño3.4 区 SST 的 3 个月滑动平均 | DJF 1950 → JJA 2026 | 919 | `SEAS YR TOTAL ANOM` |
| `sstoi.indices` | 四个 Niño 区月平均 SST 与距平（OI SST 版本） | 1982-01 → 2026-08 | 536 | `YR MON NINO1+2 ANOM NINO3 ANOM NINO4 ANOM NINO3.4 ANOM` |
| `ersst5.nino.mth.81-10.ascii` | 基于 ERSSTv5 的 Niño 指数，**气候基准期 1981–2010** | 1950-01 → 2020-12 | 852 | 同上 |
| `ersst5.nino.mth.91-20.ascii` | 基于 ERSSTv5 的 Niño 指数，**气候基准期 1991–2020** | 1950-01 → 2026-06 | 918 | 同上 |
| `Readme.index.shtml` | CPC 官方 FAQ 页面存档（HTML） | — | — | — |

**Niño 区定义**（来自 Readme）
| 区域 | 纬度 | 经度 |
|------|------|------|
| Niño 1+2 | 0–10°S | 90°W–80°W |
| Niño 3 | 5°N–5°S | 150°W–90°W |
| Niño 3.4 | 5°N–5°S | 170°W–120°W |
| Niño 4 | 5°N–5°S | 160°E–150°W |

**单位与约定**
- `TOTAL` / 各 Niño 列：区域平均 SST，单位 **°C**；
- `ANOM`：相对于基准期气候态的距平，单位 **°C**；
- ONI 的 `SEAS` 为三月重叠季（DJF、JFM、…），非日历月；该 Readme 只说明自 2021 年 1 月起 CPC 各指数统一采用 1991–2020 基准期；ONI 另按官方文档使用逐 5 年滑动更新的 30 年中心化基准期（此点非本目录文件所载，引用前请核对 CPC 官网）。
- 时间分辨率：月（ONI 为 3 个月滑动）；空间分辨率：无格点，为区域平均标量时间序列。
- 数据更新频率：每月 10 日前后（Readme）。

**注意**：`ersst5.nino.mth.81-10.ascii` 已于 **2020-12 停更**（旧基准期存档），当前业务基准期为 1991–2020（自 2021 年 1 月起启用）；日常分析应使用 `ersst5.nino.mth.91-20.ascii`。三个 Niño 文件覆盖时段重叠但基准期/生成方式不同，数值不可混用。

---

# 4. `/data1/lorraine-20419/ocean_data/ERSSTv5`

**文件**：`sst.mnmean.nc`（152 MB，单文件）
**产品**：NOAA Extended Reconstructed Sea Surface Temperature v5（in situ only），由 NOAA/ESRL PSD 重新打包为时间聚合版本。
**Conventions**：`CF-1.6, ACDD-1.3`；`processing_level: NOAA Level 4`；`product_version: Version 5`。

| 维度 | 大小 | 范围 | 分辨率 |
|------|-----|------|--------|
| `time` | 2072 | 1854-01-01 → 2026-08-01 | 月（月初打标，带 `time_bnds`） |
| `lat` | 89 | −88 → 88 °N | 2.0°（均匀） |
| `lon` | 180 | 0 → 358 °E | 2.0°（均匀） |

**变量**

| 变量 | 维度 | 单位 | 说明 |
|------|------|------|------|
| `sst` | (time, lat, lon) | degC | Monthly Means of Sea Surface Temperature；`level_desc: Surface`；`valid_range` = [−1.8, 45.0]；实际范围 [−1.8, 42.33] |
| `time_bnds` | (time, nbnds) | — | 月时间边界 |

**关键说明**
- 垂向：仅表层（无 depth 维）。
- 气候态基准：全局属性 `climatology` 声明基于 **1971–2000** SST（Xue et al., 2003）。注意这是产品内部气候态定义，文件中存的是绝对温度而非距平。
- 数据来源：现场观测 ICOADS R3.0（2015 年前）+ NCEP in situ GTS（2016 起）+ Argo SST（1999 起）；海冰 HadISST2（2015 前）/ NCEP ice（2015 后）。
- 引用：Huang et al., 2017, *J. Climate*, doi:10.1175/JCLI-D-16-0836.1。
- `data_modified: 2026-09-03`，即近实时更新的数据副本。

---

# 5. `/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies`

**产品**：GLODAPv2 全球海洋碳化学**映射气候态**（Mapped Climatologies），由 DIVA 客观分析生成。
**机构**：University of Bergen（Siv K. Lauvset 等），创建于 2016-05-12。
**引用**：Lauvset et al., *Earth Syst. Sci. Data*（GLODAPv2 mapped climatology）。
**文件数**：14 个 NetCDF，每个约 103 MB，合计 1.4 GB。

## 5.1 共同网格（所有 14 个文件完全一致）

| 维度 | 大小 | 范围 | 分辨率 |
|------|-----|------|--------|
| `lon` | 360 | 20.5 → 379.5 °E | 1.0°（均匀，**起点为 20°E**） |
| `lat` | 180 | −89.5 → 89.5 °N | 1.0°（均匀） |
| `depth_surface` | 33 | 由 `Depth` 变量给出，0 → 5500 m | 非均匀 |

`Depth`（33 个标准层，单位 m，`positive: down`）：
`0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500`

**时间**：无时间维。为气候态场，输入数据取自 **1972–2013**（各文件 `time_period` 属性一致声明：surfaces 1–33 均使用 1972–2013 年数据）。

## 5.2 文件与主变量

| 文件 | 主变量 | 单位 | 含义 |
|------|--------|------|------|
| `GLODAPv2.2016b.temperature.nc` | `temperature` | degrees celcius(sic) | 海水温度 |
| `GLODAPv2.2016b.salinity.nc` | `salinity` | 无量纲 | 实用盐度 |
| `GLODAPv2.2016b.oxygen.nc` | `oxygen` | µmol kg-1 | 溶解氧 |
| `GLODAPv2.2016b.NO3.nc` | `NO3` | µmol kg-1 | 硝酸盐 |
| `GLODAPv2.2016b.PO4.nc` | `PO4` | µmol kg-1 | 磷酸盐 |
| `GLODAPv2.2016b.silicate.nc` | `silicate` | µmol kg-1 | 硅酸盐 |
| `GLODAPv2.2016b.TCO2.nc` | `TCO2` | µmol kg-1 | 溶解无机碳（DIC） |
| `GLODAPv2.2016b.PI_TCO2.nc` | `PI_TCO2` | µmol kg-1 | 工业革命前 DIC |
| `GLODAPv2.2016b.Cant.nc` | `Cant` | µmol kg-1 | 人为碳 |
| `GLODAPv2.2016b.TAlk.nc` | `TAlk` | µmol kg-1 | 总碱度 |
| `GLODAPv2.2016b.pHtsinsitutp.nc` | `pHtsinsitutp` | 无量纲 | pH（total scale，原位温压） |
| `GLODAPv2.2016b.pHts25p0.nc` | `pHts25p0` | 无量纲 | pH（total scale，25 °C / 0 dbar） |
| `GLODAPv2.2016b.OmegaA.nc` | `OmegaA` | 无量纲 | 文石饱和度（原位温压） |
| `GLODAPv2.2016b.OmegaC.nc` | `OmegaC` | 无量纲 | 方解石饱和度（原位温压） |

## 5.3 每个文件的辅助变量（结构完全相同）

| 变量 | 维度 | 说明 |
|------|------|------|
| `<var>` | (depth_surface, lat, lon) | DIVA 客观分析后的映射场 |
| `<var>_error` | (depth_surface, lat, lon) | **映射误差**（仅内插误差，不含观测/计算不确定度） |
| `<var>_relerr` | (depth_surface, lat, lon) | 相对误差 = error / 输入数据标准差 |
| `Input_mean` | (depth_surface, lat, lon) | 输入数据经三次 Hermite 样条垂向插值后、在 1°×1° 格内的平均 |
| `Input_std` | (depth_surface, lat, lon) | 上述格内标准差 |
| `Input_N` | (depth_surface, lat, lon) | 格内数据个数 |
| `SnR` | (snr,) = 1 | 信噪比，值 = **10** |
| `CL` | (snr,) = 1 | 相关长度，值 = **7**（单位 degrees_north；纬向按 2 倍缩放，即 14°） |
| `Depth` | (depth_surface,) | 33 层深度，m |

**注意**
- `Input_mean` / `Input_std` 的 `units` 属性在**所有文件里都被硬写成 `micro-mol kg-1`**，包括 temperature、salinity、pH、Omega 等文件。这是产品的元数据错误——实际单位应与该文件主变量一致。程序化读取单位时必须以主变量的 `units` 为准。
- `temperature` 的单位字符串拼写为 `degrees celcius`（原文拼写错误），无法被 `cf_units`/`pint` 直接解析。
- 全局属性 `Comment` 字段内容为乱码（字符交错），不可用。
- 无时间维、无经度环绕元数据，与 CMIP6（0–360）或 WOA（−180–180）比对前需重排经度。

---

# 6. `/data1/lorraine-20419/ocean_data/SEANOE`

**产品**：全球海洋混合层深度（MLD）月气候态，LOPS 实验室（IFREMER / 布列塔尼大学 / CNRS / IRD）。
**方法**：对 Argo 与 NOAA NCEI World Ocean Database 的单条剖面逐条估算 MLD，再逐月做客观分析 / kriging 映射。

| 文件 | 版本 | 参考深度 | DOI | 大小 |
|------|------|---------|-----|------|
| `mld_dr003_ref10m_v2023.nc` | v2023 | **10 m** | 10.17882/91774 | 6.4 MB |
| `mld_dr003_ref5m_v2024.nc` | v2024 | **5 m** | 10.17882/98226 | 6.4 MB |
| `product_doc_v2023.pdf` | v2023 产品文档 | — | — | 19 MB |
| `product_doc_v2024.pdf` | v2024 产品文档 | — | — | 20 MB |

## 6.1 网格与时间（两个 .nc 完全一致）

| 维度 | 大小 | 范围 | 分辨率 |
|------|-----|------|--------|
| `lon` | 360 | −179.5 → 179.5 °E | 1.0°（均匀） |
| `lat` | 180 | −89.5 → 89.5 °N | 1.0°（均匀） |
| `time` | 12 | 14 → 348 | 月气候态（12 个月） |

**观测时段**：`time` 的 `climatology` 属性声明——月气候态基于 **1970-01-01 至 2021-12-31** 之间可用的剖面，时间坐标取各月中心。

**时间坐标警告**：`time` 取值为 `14, 45, 73, 104, 134, 165, 195, 226, 257, 287, 318, 348`，单位属性写作 **`unit`（非 CF 标准的 `units`）**，值为 `days since 1-1-1 0:0:0`。因此 `xarray`/CF 工具**不会自动解码**，会当作无量纲整数读入；这些数值实为各月中心的年内日序，使用前需自行映射为 1–12 月。`lat`/`lon` 同样使用 `unit` 而非 `units`。

## 6.2 变量

| 变量 | 维度 | 单位 | 说明 |
|------|------|------|------|
| `mld_dr003` | (time, lat, lon) | **m** | 混合层深度，判据：位势密度相对参考深度增加 **+0.03 kg m⁻³**。v2023 参考 10 m（有效值 10.01–1221.18 m）；v2024 参考 5 m（有效值 5.00–1236.94 m） |
| `n_obs` | (time, lat, lon) | 计数 | 每格每月的现场观测（剖面）数量；v2023 最大 27760，v2024 最大 27792 |
| `mask` | (lat, lon) | 0/1 | 陆海掩膜：0=陆，1=海。由 etopo2 双线性插值到 0.1° 生成 |

**数据来源**：Argo 计划 + NOAA NCEI World Ocean Database，剖面类型涵盖 CTD、APB（海洋哺乳动物标记）、OSD（采水瓶）、MRB（锚系浮标，含 TAO 阵列）、GLD（滑翔机）、UOR（波动式拖曳记录仪）、DRB（漂流浮标，含 Ice Tethered Profiler）。

**注意**
- 所有坐标与数据变量的单位属性均写作 `unit`（非 CF 标准的 `units`），CF 解析工具不会自动识别。
- v2023 与 v2024 的参考深度不同（10 m vs 5 m），MLD 数值下限随之不同，**两版不可直接互换或做差**。
- v2024 文件的 `Timestamp` 属性为 `Fri Jul 14 15:15:48 2023`，早于版本号所示年份，属产品元数据笔误。
- 无垂向维度（MLD 本身即深度标量场）。

---

# 7. `/data1/lorraine-20419/ocean_data/WOA23`

**产品**：World Ocean Atlas 2023（NOAA/NCEI），1.00° 客观分析气候态。
**Conventions**：`CF-1.6` + `NCEI_NetCDF_Grid_Template_v1.0`；`crs` 变量声明 EPSG:4326。
**目录结构**：`climatology/`（6 个变量子目录 × 17 文件 = 102 个 .nc，3.3 GB） + `masks/`（4 个 ASCII 掩膜，171 MB）。

## 7.1 文件命名与时间片编号

`woa23_<period>_<var><NN>_01.nc`
- `<period>`：`decav71A0` = 1971–2000 十年平均；`all` = 全部可用年份
- `<var>`：`t`=温度, `s`=盐度, `o`=溶解氧, `n`=硝酸盐, `p`=磷酸盐, `i`=硅酸盐
- `<NN>`：`00`=年平均；`01`–`12`=1–12 月；`13`=冬(1-3月), `14`=春(4-6月), `15`=夏(7-9月), `16`=秋(10-12月)
- 末尾 `_01` = 1.00° 网格

## 7.2 六个变量子目录

| 子目录 | 变量前缀 | 主变量单位 | 标准名 | 气候期 | 文件数 |
|--------|---------|-----------|--------|--------|-------|
| `temperature_1971-2000` | `t` | degrees_celsius | sea_water_temperature | 1971-01-01 起 P30Y（1971–2000） | 17 |
| `salinity_1971-2000` | `s` | 1（实用盐度） | sea_water_practical_salinity | 1971-01-01 起 P30Y（1971–2000） | 17 |
| `oxygen_all_time_periods` | `o` | **µmol kg⁻¹** | moles_of_oxygen_per_unit_mass_in_sea_water | 1965-01-01 起 P58Y（1965–2022） | 17 |
| `nitrate_all_time_periods` | `n` | **µmol kg⁻¹** | moles_of_nitrate_per_unit_mass_in_sea_water | 1965–2022 | 17 |
| `phosphate_all_time_periods` | `p` | **µmol kg⁻¹** | moles_of_phosphate_per_unit_mass_in_sea_water | 1965–2022 | 17 |
| `silicate_all_time_periods` | `i` | **µmol kg⁻¹** | moles_of_silicate_per_unit_mass_in_sea_water | 1965–2022 | 17 |

**温盐与营养盐/氧的气候期不同**（1971–2000 vs 1965–2022），联合使用（如计算 AOU、密度加权营养盐）时须说明这一时段错配。

## 7.3 水平网格（所有文件一致）

| 维度 | 大小 | 范围 | 分辨率 |
|------|-----|------|--------|
| `lat` | 180 | −89.5 → 89.5 °N | 1.00°（均匀，带 `lat_bnds`） |
| `lon` | 360 | **−179.5 → 179.5** °E | 1.00°（均匀，带 `lon_bnds`） |

## 7.4 垂向网格（随变量与时间片变化——最易踩坑处）

| 文件类型 | 层数 | 最深 | 说明 |
|---------|-----|------|------|
| 温 / 盐 / 氧 / 营养盐 **年平均**（`00`） | **102** | 5500 m | 完整标准层 |
| 温 / 盐 **季节**（`13`–`16`） | **102** | 5500 m | 完整标准层 |
| 温 / 盐 / 氧 **月**（`01`–`12`） | **57** | 1500 m | 仅上层 |
| 氧 **季节**（`13`–`16`） | **57** | 1500 m | 仅上层 |
| 硝酸盐 / 磷酸盐 / 硅酸盐 **月 + 季节**（`01`–`16`） | **43** | 800 m | 仅上层 |

深度层（单位 m，`positive: down`，43 / 57 / 102 层为逐级前缀嵌套）：
- **0–100 m**：0, 5, 10, …, 100（Δ=5 m，共 21 层）
- **100–500 m**：125, 150, …, 500（Δ=25 m）
- **500–2000 m**：550, 600, …, 2000（Δ=50 m）
- **2000–5500 m**：2100, 2200, …, 5500（Δ=100 m）
- 43 层 = 至 800 m；57 层 = 至 1500 m；102 层 = 至 5500 m。

## 7.5 时间轴

`time` 长度均为 1，单位 `months since 1971-01-01 00:00:00`（**所有文件统一使用该 epoch，包括 1965–2022 的营养盐文件**），附 `climatology_bounds` 给出气候期起止。

> **警告**：该单位（`months since ...`）**不是 CF 可解析的时间单位**，`xarray.open_dataset()` 默认会报错。必须使用 `decode_times=False` 打开，并自行按 `climatology_bounds` 与文件名中的 `<NN>` 解释时间片。

## 7.6 每个文件的变量族（`<v>` 为变量前缀 t/s/o/n/p/i）

| 变量 | 单位 | 含义 |
|------|------|------|
| `<v>_an` | 同主变量 | **客观分析平均场**（最常用） |
| `<v>_mn` | 同主变量 | 统计平均（各格内所有有效插值观测的均值） |
| `<v>_dd` | 1 | 各格各标准层的观测个数 |
| `<v>_sd` | 同主变量 | 统计均值的标准差 |
| `<v>_se` | 同主变量 | 统计均值的标准误 |
| `<v>_oa` | 同主变量 | 统计均值 − 客观分析均值（差值） |
| `<v>_gp` | 1 | 最小影响半径内含统计均值的格点数 |
| `<v>_sea` | 同主变量 | 客观分析场的标准误（Levitus et al., 2012 GRL 定义） |
| `<v>_sdo` | 同主变量 | 客观分析标准差场（**仅 t / s / o 有；n / p / i 无**） |
| `<v>_ma` | 同主变量 | 月/季文件专有（**年平均文件 `00` 无此变量**） |
| `crs` | — | 坐标参考系：latitude_longitude, EPSG:4326, WGS84 椭球 |
| `lat_bnds` `lon_bnds` `depth_bnds` `climatology_bounds` | — | 各维边界 |

变量个数：年平均 t/s/o = 9 个；年平均 n/p/i = 8 个；月/季 t/s = 10 个；月/季 o = 9 个；月/季 n/p/i = 9 个。

## 7.7 `masks/` 子目录（ASCII，CSV 格式）

| 文件 | 大小 | 行数 | 网格 | 首行注释 | 列 |
|------|------|-----|------|---------|-----|
| `landsea_01.msk` | 1.2 MB | 64802 | 1.00° | `#Land/Sea mask: 1.00 grid : WOA13` | `Latitude, Longitude, Bottom_Standard_Level` |
| `landsea_04.msk` | 18 MB | 1036802 | 0.25° | `#Land/Sea mask: 0.25 grid : WOA13` | `Latitude, Longitude, Bottom_Standard_level` |
| `basinmask_01.msk` | 8.9 MB | 41090 | 1.00° | `#Basin number for WOA23 climatology objective analysis : 1.00 degree grid` | `Latitude, Longitude, Basin_0m … Basin_5500m`（102 个深度列） |
| `basinmask_04.msk` | 143 MB | 677826 | 0.25° | 同上，0.25° | 同上（102 个深度列） |

**说明**
- 陆海掩膜的 `Bottom_Standard_Level` 为该格点的**最深标准层序号**（1 = 仅表层即陆地/极浅，最大 102 = 5500 m），非 0/1 二值。
- 海盆掩膜按 102 个标准深度逐层给出海盆编号（如 10 = 南大洋，11 = 北冰洋 等 WOA 海盆编码）；超过该点海底深度的层为**空字段**（连续逗号）。
- 掩膜经度约定为 **−180–180**，与 climatology 的 NetCDF 一致；但 `landsea_01.msk` 的经度按 0.5 → 359.5 再折回负值排列（末行为 −0.5），读取时按 `Latitude, Longitude` 两列显式重建索引，不要假设行序。
- `basinmask_*` 只含**有海洋格点的行**（1° 版 41090 行 < 180×360=64800），是稀疏表；`landsea_*` 为完整格点表（64802 = 64800 + 2 行表头/注释）。
- 掩膜标注来源为 **WOA13**（沿用），并非 WOA23 重新生成。

---

# 8. 跨数据集使用清单（Checklist）

打开 / 对齐这些数据时，按顺序确认：

1. **经度约定**：CMIP6 & ERSSTv5 = 0–360；WOA23 & SEANOE = −180–180；GLODAP = 20.5–379.5。
2. **深度单位**：CESM2 = cm，其余 = m。CESM2 需 `lev / 100`。
3. **网格类型**：CESM2 与 MPI-ESM1-2-LR 是曲线网格（二维 lat/lon + 索引维），不能按 `lat`/`lon` 直接 `.sel()`；GFDL-ESM4 已在 `gr` 上重网格为规则 1°，可直接切片。
4. **时间解码**：WOA23 必须 `decode_times=False`（`months since`）；SEANOE 的 `time` 无单位；CMIP6 三种日历各异。
5. **变量缺口**：ssp245 全模式缺 `tos`/`mlotst`；ssp245/GFDL-ESM4 缺全部 `Ofx` 固定场（借用 historical）。
6. **成员一致性**：`historical/CESM2/phycos`(r6) 与 `historical/MPI/ph`(r10) 与同目录其他变量不是同一成员。
7. **单位属性可信度**：GLODAP 的 `Input_mean`/`Input_std` 单位属性一律错标为 `micro-mol kg-1`；GLODAP 温度写作 `degrees celcius`；SEANOE 用 `unit` 而非 `units`。
8. **气候期错配**：WOA23 温盐为 1971–2000，营养盐/氧为 1965–2022；GLODAP 为 1972–2013；ERSSTv5 内部气候态为 1971–2000；CPC 各指数自 2021-01 起用 1991–2020（ONI 另用滚动 30 年中心化基准期）。
9. **摩尔浓度基准不同**：CMIP6 用 **mol m⁻³**（体积基），WOA23 与 GLODAP 用 **µmol kg⁻¹**（质量基）。互转需要原位密度：`mol m⁻³ = µmol kg⁻¹ × ρ / 1000`。
