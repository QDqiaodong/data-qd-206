# 测光发布件与变暗判据契约（build_photometric_release / assess_fade）

本契约在既有 `process_exposures` 几何关联与锚点规则之上，新增两个只依赖标准库的离线 API。所有星等、零点、差值、相等与大小关系均按输入十进制的精确有理值（`Fraction`）计算，不用浮点误差或四舍五入放行；对外数值一律是十进制字符串或 `null`。

```python
from star_crossmatch import build_photometric_release, assess_fade
```

## 1. 输入 `build_photometric_release(data)`

`data` 恰含四项，字段严格、类型严格，任何非法输入抛 `ValueError`；不修改调用者数据，不交付部分结果。

| 字段 | 规模 | 行结构 |
| --- | --- | --- |
| `photometry` | 1..32 次曝光、1..256 个参考源、全批 1..2000 条检测 | 现有 `process_exposures` 输入；新增规模上限先校验，原曝光字段/类型/坐标/ID 校验照常 |
| `epochs` | 每次曝光恰一行 | `{exposure, mjd}`，MJD 为非负十进制天数；不同曝光允许同刻 |
| `transfers` | 0..128 行 | `{id, from_exposure, to_exposure, delta}`，`id` 批内唯一，两端为不同的本批曝光 |
| `limits` | 0..5000 行 | `{exposure, reference, limiting_instrumental_magnitude}`，每个（曝光，参考源）组合最多一行 |

ID 与十进制字符串语法沿用仓库统一约定（`[A-Za-z0-9_:-]{1,64}`；最多 9 位整数、6 位小数的有界十进制）。

`delta` 是有符号星等差，方向固定为 **后一字段的零点减前一字段的零点**，即约束

```
zero[to_exposure] - zero[from_exposure] = delta
```

与时间顺序无关。同一对帧之间允许多条独立记录，它们是不同的排错单位，**不预先合并**。数据生产方已完成同滤光带、稳定源与记录筛查；本 API 把差值视为精确约束，不重新从图像或检测估算。

绝对基准（锚点）是本批硬约束，不在故障模型中：仅当某帧沿用原锚点规则得到有效零点时，该帧才是绝对定标基准，基准不可被“排除”。由转移得到的零点不反过来改变锚点，也不用于确认检测身份。零点无效的帧没有基准。

## 2. 两类测光证据的核验

### 2.1 可信实测（检测身份）

只使用曝光内的**几何最优方案**（先最大化配对数、再最小化整数角距代价之和，与 `process_exposures` 同一套最小费用流/零费用残差环流刻画），完全不做星等破平局：

* 对某参考源，只有在**所有**几何最优方案中都配到**同一条**检测（对应边 `flow=1` 且不存在令其为 0 的最优方案）时，才把该检测的**原始仪器星等**签发为该（曝光，参考源）的可信实测；
* 不能用有效零点下的星等残差破平局来签发实测。`process_exposures` 能靠星等选定的配对，若几何上并列，发布件中仍不是可信实测。

### 2.2 灵敏度证明（未检出证据）

`limits` 是上游在未遮挡位置做的灵敏度证明：若源在该帧不比 `limiting_instrumental_magnitude` 暗（仪器星等越大越暗），原则上会留下提取记录。

* 只要该参考源出现在该帧**任一检测的原始几何候选**里——包括该候选因更优分配从未成为最终配对的情况——该证明就被阻止，不能作为未检出证据；发布件保留这些阻止检测的 ID；
* 没有 `limits` 行意味着没有这份证明，不能从“曝光存在、别的源有深度、未选中配对”推断；
* 证明本身不说明是否实际检出。只有通过上述核验的清洁未检出，才可推出该帧该源的星等**严格大于**极限值。

为使早于定标的拒绝可追溯，发布件另含 `geometric_identities`：对“有原始几何候选但未锁定可信实测”的（曝光，参考源）组合，仅保留候选检测 ID 列表；不含坐标、角距或候选边表。

## 3. 零点能发布到什么程度

### 3.1 无向连通分量与假设枚举

先按**全部**转移记录忽略方向划分无向连通分量；没有任何边的曝光自成一个分量。每个分量独立保留各自的可行假设：

1. 若全部转移与本分量的绝对基准可以同时成立，只保留 `excluded_transfer = null` 的“不删记录”假设一种；
2. 否则尝试**恰好排除一条**转移记录，保留所有能使本分量**全部剩余约束**（含不可删除的绝对基准）成立的单条排除结果；
3. 不能用 ID 只挑一个能通的假设；需要排除两条及以上的分量没有可行假设（`status="unsolvable"`，假设列表为空），不得当成已修复。

各分量的故障预算与假设互不影响；任一分量无解不抹掉别组的成果。

### 3.2 假设内的确定量

一个假设内仍互通的各帧有确定的零点差；仅当互通子组中含绝对基准时，组内帧才有确定绝对值。没有基准的子组可以整体平移，不默认零点为 0。

每帧仅在本分量**每个**可行假设都给出**同一个**确定绝对值时才发布该值（`zero_point_status="published"`），否则 `zero_point` 为 `null`，状态区分：

| 状态 | 含义 |
| --- | --- |
| `published` | 全部可行假设给出同一确定绝对值 |
| `indeterminate` | 分量有可行假设，但该帧绝对值不存在或假设间不一致 |
| `unanchored` | 该帧所在互通结构在任何假设下都没有绝对基准（可整体平移） |
| `unsolvable` | 所属分量没有可行假设（含该帧本身是锚点的情形） |

## 4. 发布件格式（版本化 JSON）

`build_photometric_release` 返回 `PhotometricRelease`（`dict` 子类），其 dict 本身即载荷，`format = "star-crossmatch-photometric-release/1"`，可直接 `json.dumps`；也可用 `.to_json()`。顶层字段：

* `format`：版本标识；未知格式以 `ValueError` 拒绝。
* `references`：参考源 ID 排序列表（只保留身份）。
* `exposures`：按曝光 ID 排序，每项 `{exposure, mjd, is_anchor, zero_point_status, zero_point}`。
* `transfers`：按转移 ID 排序，每项 `{id, from_exposure, to_exposure, delta, component}`。
* `components`：按分量成员身份排序，每项 `{component, exposures, status, hypotheses}`；`hypotheses` 每项 `{hypothesis, excluded_transfer}`，假设 ID 形如 `component-1#h1`，排序规则为 `null` 排除优先、其余按被排除转移 ID 字典序。
* `measurements`：可信实测，每项 `{exposure, reference, observation, instrumental_magnitude}`，按（曝光，参考源）排序。
* `sensitivity_proofs`：每项 `{exposure, reference, limiting_instrumental_magnitude, usable, blocking_detections}`，按（曝光，参考源）排序。
* `geometric_identities`：每项 `{exposure, reference, candidate_detections}`，按（曝光，参考源）排序。

发布件**不**包含原始曝光行、坐标、匹配候选边表或整份 `photometry`；被阻止使用未检出证明的组合只保留检测 ID。重排任何输入行都不改变发布件（全部输出有确定排序）。

装载（含跨进程的纯 dict）时会复算每个假设的全部约束并校验状态/身份一致性，篡改出的不可行假设或错误分量归属以 `ValueError` 拒绝。

## 5. `assess_fade(release, query)`

查询为 `{reference, bright_exposure, later_exposure, minimum_fade}`，`minimum_fade` 为正十进制星等。发布件可以是 `PhotometricRelease`，也可以是经 JSON 往返的纯 dict；新进程无需原始输入、图像或任何缓存。

未知身份、坏数据、不支持的发布件格式、非正幅度或前帧 MJD 不严格小于后帧（含同刻）抛 `ValueError`。物理上的定标冲突与证据不足是正常业务结果。

### 5.1 判定顺序与结果状态

1. 观测证据门（早于定标，状态互不混淆，结论带曝光、参考源，适用时带检测 ID）：
   * `no_trusted_measurement`：前帧没有覆盖全部几何最优方案的可信实测；`detections` 给出 `geometric_identities` 中适用的候选检测 ID（可能为空）。
   * `no_sensitivity_proof`：后帧没有该源的灵敏度证明。
   * `sensitivity_proof_blocked`：后帧存在该源的原始几何候选，证明不可用；`detections` 明确列出阻止它的检测 ID。
2. 定标门：查询所用任一分量 `unsolvable`，返回 `calibration_unsolvable`。
3. 联合假设：
   * 两帧在同一原分量时，必须在**同一个**排除假设下作差，不能把两帧各自的独立上下界拼差；
   * 两帧在不同原分量时，两组假设可以任意组合；
   * 任一联合假设中两帧不互通且至少一方无确定绝对值（就该假设而言），零点差不受约束，返回 `not_comparable`，并给出一个达到该情形的假设证据；两帧虽不互通但该假设下都有确定绝对值时仍可比较。
4. 下界：对所有可行且可比的联合假设，计算

   ```
   bound = 后帧仪器极限 − 前帧实测仪器星等
           + 后帧零点 − 前帧零点
   ```

   返回其中的最小值作为 `guaranteed_fade_lower_bound`（十进制字符串）。`bound >= minimum_fade` 时 `status="fade_confirmed"`，否则 `status="insufficient_evidence"`（不表示源没有变暗）。

每个可比结论都带 `hypothesis_evidence`：涉及的分量及该假设排除的转移记录（未排除为 `null`）。并列时的文档固定规则：联合假设按分量 ID、再按发布件中的假设固定次序枚举，取首个达到最小下界者（即字典序最小的分量/排除组合）。

## 6. 离线交接

```python
import json
release = build_photometric_release(data)
blob = json.dumps(release, ensure_ascii=False, sort_keys=True)
# 交接方:新进程,只有 blob
verdict = assess_fade(json.loads(blob), query)
```

`tests/test_release.py` 用 `subprocess` 以 `python -I` 启动独立子进程，只经 stdin 传入发布件 JSON 与查询，验证子进程结论与原进程逐项一致。`examples/run_release.py` 是可离线运行的最小演示。
