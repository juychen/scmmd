# Region-level subset loom（`11subsetloom_meth.py`）

## 1. 原有 loom 提取算法（本目录的做法）

| 脚本 | 输入 | 粒度 | 关键点 |
|---|---|---|---|
| `11subsetloom.py` | `stg1028/{CUMS,CURES,CSRES,CSDS}_4VN/*.h5ad` | `region + celltype.L2` | 逐个数据集按 `status` 过滤，拼 `celltype.L2.region` 名字，`anndata.concat(join="inner")` 后 `loompy.create` |
| `11subsetloom_new.py` | 单个 region 的 h5ad（如 `/data3/junyi/scvi/PFC.h5ad`） | `region + celltype` | 支持 `--layer`，单文件 |
| `11subsetloom_3v3.py` | 多个 3v3 降采样 h5ad | `region + celltype` | 用 obs 名去重共享 control 细胞，`lp.create` |

共同算法（也就是新脚本继承的部分）：

1. `anndata.read_h5ad(..., backed="r")` 只读 obs，不载入矩阵；
2. 按 obs 的 `region` / `status` / `sex` / `celltype` 生成选择掩码；
3. 对每个输出单元，从各输入里取子集 → 拼接（`join="inner"`）；
4. 组装 `row_attrs = {"Gene": var_names}`，
   `col_attrs = {"CellID", "nGene", "nUMI"}`；
5. `loompy.create(path, X.T, row_attrs, col_attrs)`（行=基因，列=细胞）。

新脚本保留 1–4 的思路，但改成 **一个脑区一个 loom**（而不是一个 celltype 一个 loom），
并且换了 3 处实现，见第 4 节。

## 2. 分组设计（默认 `GROUP_SPECS`）

| Group | 数据集 | `status` | `sex` | 脑区（obs 原值） |
|---|---|---|---|---|
| `CUMS_M` | CUMS_4VN | `SUS` | M | AMY PFC TH **HPF** |
| `CUMS_F` | CUMS_4VN | `SUS` | F | AMY PFC TH **HPF** |
| `CSDS_M` | CSDS_4VN | `CSDS` | M | AMY PFC TH **HPF** STR HY MB |
| `CON_M` | CUMS_4VN | `CON` | M | AMY PFC TH **HPF** STR HY MB |
| `CON_F` | CUMS_4VN | `CON` | F | AMY PFC TH **HPF** |

> 数据里海马叫 `HPF`，输出默认重命名为 `HIP`（`--rename HPF=HIP`），
> 与 ATAC 那边的命名保持一致；obs 原值仍记在 `RegionSource`。

输出（一个脑区一个 loom）：

```
AMY.loom / PFC.loom / TH.loom / HIP.loom   ← CUMS_M + CUMS_F + CSDS_M + CON_M + CON_F
STR.loom / HY.loom  / MB.loom              ← CSDS_M + CON_M
```

细胞数（实测，未降采样）：

| region | CUMS_M | CUMS_F | CSDS_M | CON_M | CON_F | total |
|---|---|---|---|---|---|---|
| AMY | 42,933 | 41,650 | 24,564 | 44,028 | 33,107 | 186,282 |
| PFC | 33,857 | 48,709 | 27,110 | 37,463 | 38,034 | 185,173 |
| TH  | 32,126 | 26,056 | 23,234 | 27,243 | 27,409 | 136,068 |
| HIP | 38,606 | 36,239 | 32,036 | 41,373 | 26,765 | 175,019 |
| STR | – | – | 27,508 | 33,104 | – | 60,612 |
| HY  | – | – | 24,129 | 37,449 | – | 61,578 |
| MB  | – | – | 34,839 | 40,640 | – | 75,479 |

**control 去重**：`CSDS_4VN` 里的 289,670 个 `CON` 细胞与 `CUMS_4VN` 的 `CON` 男性细胞
完全重合（已核对 barcode）。脚本按 barcode 全局去重，所以同一个 control 细胞不会进两次。

## 3. loom 内容

`row_attrs`：`Gene`

`col_attrs`：

| 列 | 说明 |
|---|---|
| `CellID` | h5ad 的 obs 名（含 sample 后缀，全局唯一） |
| `Group` | `CUMS_M` / `CUMS_F` / `CSDS_M` / `CON_M` / `CON_F` |
| `Condition` | `CUMS` / `CSDS` / `CON` |
| `Sex` | `M` / `F` |
| `Dataset` | `CUMS` / `CSDS` |
| `Region` | 输出脑区名（`HIP`） |
| `RegionSource` | obs 原值（`HPF`） |
| `Sample`, `Donor` | 原 obs 列 |
| `celltype.L1`, `celltype.L2` | 原 obs 列 |
| `nGene`, `nUMI` | 由写入的矩阵计算 |

矩阵默认取 `layers["counts"]`（原始 UMI，和 `sn0827L1/*.loom` 一致）；`--layer X` 取归一化矩阵。

## 4. 与旧脚本不同的 3 处实现

1. **读 obs / 矩阵不再用 `anndata` backed**
   CUMS 的 obs 是 1.1 M × 69 列，`read_h5ad(backed="r")` 要几分钟、几个 GB。
   现在用 `h5py` 直接读：分类型列只读 `codes`（1 字节/细胞），需要时再解码；
   矩阵按 CSR 的 `indptr` 做行切片（`CsrSource`）。
2. **拼接不再用 `scipy.sparse.vstack`**
   改成 `vstack_csr()`：一次性预分配 `data/indices/indptr` 再逐块 memcpy。
3. **写 loom 不再用 `loompy.create`**
   `loompy.create` 写稀疏矩阵时把列窗口 densify 到 **(64,64) 的 gzip chunk** 里，
   区域级矩阵（16k × 186k、4.5 亿非零）会产生上百万个小 chunk，实测卡住十几分钟不出文件。
   现在 `create_loom_fast()` 用 h5py 直接写 `(1024,1024)` chunk 的 `/matrix`，
   仍是标准 Loom 3.0 文件，`loompy.connect` / pySCENIC / SCopeLoomR 读取无差别
   （已用 `loompy.connect(validate=True)` 和逐细胞比对验证）。

## 5. 用法

```bash
# 只看选择结果（不写文件，几秒钟）
python 11subsetloom_meth.py --dry-run

# 生成全部脑区
bash 11subset_loom_meth.sh

# 只做某几个脑区；每组每脑区最多 2 万细胞
bash 11subset_loom_meth.sh --regions AMY,STR --max_cells_per_group 20000

# 重建（覆盖已有 loom）
bash 11subset_loom_meth.sh --overwrite

# 只检查已有 loom 并重写 manifest
python 11subsetloom_meth.py --manifest-only

# 换输入 / 换分组
python 11subsetloom_meth.py \
    --input CUMS=/path/CUMS.h5ad --input CSDS=/path/CSDS.h5ad \
    --status_column status --sex_column sex --region_column region \
    --output /data1st1/junyi/output/meth_region_loom
```

`--config config.json` 可以一次给出 `{"inputs": {...}, "groups": [...], "rename": {...}}`，
`groups` 的结构和文件里的 `GROUP_SPECS` 相同（`dataset/status/label/sex/regions`）。

输出目录里会写：

* `*.loom`：每个脑区一个；
* `loom_manifest.csv`：region / status / n_cells / n_genes / size / groups；
* `run_config.json`：本次运行的输入、分组、参数。

## 6. 运行开销（实测，384 核 / 1.5 T 内存机器）

* 首次读入两个 h5ad 的矩阵缓存：CUMS 20.3 GB ≈ 4.5 min，CSDS 9.2 GB ≈ 1.3 min（NFS）；
* 每个脑区切片 ≈ 5–35 s，写 loom ≈ 35–100 s；
* 7 个脑区全量一次跑完约 25 min（含两次缓存加载），峰值内存 ~50 GB；
* 输出：880,211 个细胞 × 16,428 基因，7 个 loom 共 **4.3 GB**（gzip）。

## 7. 校验

`loompy.connect(..., validate=True)` 通过；每个 loom 的
`Group` 细胞数与预期一致、`CellID` 无重复、基因顺序与 h5ad 的 `var` 完全一致；
每个 group 随机抽 1–2 个细胞，把 loom 里的整列向量与源 h5ad `layers["counts"]` 逐值比对，完全相同。

```bash
# 重新检查已有 loom 并刷新 manifest（不重建）
python 11subsetloom_meth.py --manifest-only
```
