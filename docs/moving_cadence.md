# 三帧移动巡天归约契约（reduce_moving_cadence）

新增 `reduce_moving_cadence(data)`：一次接收同一小天区的恰好三帧检测、实际 MJD、统一平面坐标和移动速度/容差，返回：

1. 精确的移动轨迹候选、冲突分量、全批最佳值和逐条检测审计；
2. 每个分量的无约束最佳见证，以及硬要使用/不用某检测的约束最佳见证；
3. 按移动状态实际剔除或限制锚点后的静态测光结果，其候选、代价、校准和星等破平局规则与 `process_exposures` 一致。

只依赖 Python 标准库。所有输入数值仍使用仓内既有十进制字符串语法（可选负号、最多 9 位整数和 6 位小数）；错误输入抛 `ValueError`，不修改调用方对象。

## 调用

```python
from star_crossmatch import reduce_moving_cadence

result = reduce_moving_cadence(data)
static = result["static_photometry"]   # {"exposures": [...]}，结构同 process_exposures
audit = result["detections"]
```

可离线运行：

```sh
python3 -I examples/run_moving_cadence.py
```

## 输入

`data` 必须恰好包含以下字段，不接受别名或额外字段：

| 字段 | 约束 |
| --- | --- |
| `photometry` | 完整的 `process_exposures` 输入对象；本接口额外要求恰好 3 帧、每帧 1..120 条检测、参考源 1..256 个。 |
| `epochs` | 恰好 3 行，每行恰为 `{"exposure", "mjd"}`；曝光必须与 `photometry.exposures` 一一对应。MJD 为互不相等的非负十进制天数。 |
| `positions` | 每条测光检测恰好一行 `{"exposure", "observation", "x", "y"}`，不得重复、缺失或指向未知身份。`x/y` 为角秒平面坐标。 |
| `min_speed` | 正十进制字符串，角秒/天。 |
| `max_speed` | 正十进制字符串，且 `min_speed <= max_speed`。 |
| `middle_tolerance` | 非负十进制字符串，角秒。 |

检测完整身份是二元组 `(exposure, observation)`。不同曝光中允许出现相同的 `observation` 字符串。早、中、晚帧按 MJD 决定，与 `epochs` 或曝光数组顺序无关。

## 轨迹和精确代价

对每条早帧检测 `E`、中帧检测 `M`、晚帧检测 `L`：

- 令实际间隔为 `dt = mjd_L - mjd_E`；
- 平面速度长度为 `sqrt((x_L-x_E)^2 + (y_L-y_E)^2) / dt`，必须落在闭区间 `[min_speed, max_speed]`；
- 按实际 MJD 线性插值出中帧预测位置；
- `M` 到预测点的欧氏距离必须 `<= middle_tolerance`；
- 候选代价为该距离的平方。

速度边界、容差边界、入围比较和总代价均使用 `fractions.Fraction`/整数平方比较，不使用取整速度或浮点 epsilon。代价可能不是有限小数，因此输出为约分有理数字符串：整数输出 `"1"`，否则输出 `"1/100"`、`"1/200"` 这种形式。

## 冲突分量和规模上限

两条候选只要共享任一检测即相连；连通闭包构成冲突分量。分量内参与过任意候选的不同检测最多 15 条；从未进入候选的检测不计。超过 15 条时整批抛 `ValueError`，不会把过密分量静默当成无轨迹。

不同分量相互独立。每个分量先最大化轨迹条数，再最小化总代价。全批最佳 `(轨迹数, 总代价)` 是各分量最佳值之和；接口不枚举许多独立分量的笛卡尔积。

## 输出结构

下面给出完整规范字段。顶层同时提供 `epochs`，其值与 `frames` 相同（按早、中、晚排序并带 `role`）。

### 顶层

```json
{
  "schema": "star_crossmatch.moving_cadence/v1",
  "frames": [
    {"role": "early", "exposure": "E", "mjd": "60000"},
    {"role": "middle", "exposure": "M", "mjd": "60001"},
    {"role": "late", "exposure": "L", "mjd": "60002"}
  ],
  "parameters": {
    "min_speed": "0.75",
    "max_speed": "1.25",
    "middle_tolerance": "0.2"
  },
  "global_best": {"trajectory_count": 1, "total_cost": "1/100"},
  "components": ["..."],
  "candidates": ["..."],
  "detections": ["..."],
  "static_photometry": {"exposures": ["..."]},
  "static_handoff": {
    "removed_detections": ["..."],
    "anchor_blocked_detections": ["..."],
    "retained_detections": ["..."]
  }
}
```

完整静态结果在 `static_photometry`（同值别名 `photometry`），其 `exposures` 结构与 `process_exposures()["exposures"]` 完全相同；顶层 `exposures` 也是同一曝光列表。

### 身份对象

所有审计和见证中的检测身份均为：

```json
{
  "exposure": "M",
  "observation": "B",
  "identity": ["M", "B"]
}
```

`identity` 就是 JSON 数组形式的完整身份二元组。

### 分量和候选

`components[]` 含：

- `id` / `component_id`：稳定字符串，如 `C1`；
- `component_index`：从 1 开始的序号；
- `detections`：分量内参与候选的检测身份；
- `candidate_count`、`candidates`；
- `unconstrained_best`：该分量无约束最佳方案；`best` 是同值别名。

候选记录含 `id`（如 `C1-T1`）、`early`、`middle`、`late`、`cost` 和精确的 `speed_squared`。

方案对象形如：

```json
{
  "feasible": true,
  "trajectory_count": 1,
  "total_cost": "1/100",
  "witness": [
    [
      {"exposure": "E", "observation": "A"},
      {"exposure": "M", "observation": "B"},
      {"exposure": "L", "observation": "C"}
    ]
  ],
  "identity_witness": [[["E", "A"], ["M", "B"], ["L", "C"]]]
}
```

`trajectory_count`/`total_cost` 也分别提供短别名 `count`/`cost`。`witness` 中每条轨迹按早、中、晚排列；同一见证不会复用检测。`identity_witness` 是同一见证的紧凑二元组形式。并列时选择按完整身份排序后字典序最小的见证，因此不随输入行序变化。该见证是一份可复核样本，不表示其他并列最佳方案不存在。

### 检测审计

`detections`（同值别名 `audit`）按 `(exposure, observation)` 排序。参与过候选的检测含：

- `component` / `component_id`；
- `state`：
  - `"always"`：每个无约束全批最佳方案都使用；
  - `"partial"`：仅部分最佳方案使用；
  - `"never"`：参与过候选，但没有任何最佳方案使用；
- `status`：分别为 `"always_used"`、`"partially_used"`、`"never_used"`；
- `counterfactuals.require_used`：硬要使用该检测时，其分量能达到的最佳方案；
- `counterfactuals.forbid_used`：硬要不用该检测时的最佳方案。

`require_used`（别名 `if_used`）表示硬要使用该检测，`forbid_used`（别名 `if_unused`）表示硬要不用该检测。约束方案可以比无约束最佳更差；条件无解时对应值为 `null`。

从未进入候选的检测：

- `state` 为 `"not_candidate"`，`status` 为 `"not_in_candidate"`；
- `component`、`component_id` 和 `counterfactuals` 为 `null`；
- 不虚构分量或见证。

### 静态交接状态

检测记录还给出实际交接动作：

- `state="always"`：`removed_from_static=true`、`static_action="remove"`；
- `state="partial"`：`anchor_blocked=true`、`static_action="retain_but_block_anchor"`；
- 其他保留检测：`static_action="retain"`。

静态处理严格按以下顺序重新计算：

1. 从 `photometry` 中删除所有 `always` 检测，而不是事后给旧结果加删除标记；
2. 对保留检测重算原球面候选、半径和整数毫角秒代价；
3. 重算一对一几何最优和强制配对；
4. 只有未被移动状态禁止的几何强制 calibrator 配对才能进入零点锚点集合；
5. 按原锚点数量、中位数、容差支持和多数规则确定零点；
6. 按原规则计算校正星等，并在零点有效时做几何最优内的星等残差破平局。

`partial` 检测只是其自身配对不得作为锚点；限制不转移到参考星。同一参考星在其他曝光中，或在同一曝光被另一条未受限检测强制配到时，仍可正常成为锚点。

剔除后没有检测的曝光仍出现在 `static_photometry.exposures` 中：`detections=[]`，校准为 `insufficient_anchors`、`zero_point=null`、`anchors=[]`。

当不存在任何 `always` 或 `partial` 限制时，静态部分直接由 `process_exposures(data["photometry"])` 给出，因此与原输出完全一致。

## 确定性和错误

- 所有行按曝光 ID、检测 ID、MJD 或完整身份的规范排序处理；重排 `epochs`、`positions`、参考源、曝光或检测行不改变 JSON 结果。
- 分量 ID 按分量成员的完整身份元组排序后分配。
- 严格字段、范围、重复身份、缺失身份、非法十进制数和分量超限均抛 `ValueError`。
- 成功结果只含 `dict`、`list`、字符串、整数、布尔和 `null`，可直接 `json.dumps`。
