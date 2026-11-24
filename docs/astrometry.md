# 二维仿射天体测量解算契约（solve_astrometry）

`star_crossmatch.solve_astrometry(data)` 为一次曝光建立探测器 `(x,y)` 到以
`tangent_point` 为中心的 gnomonic 切平面 `(xi,eta)` 的六参数二维仿射解，再把
科学检测反投影到 RA/Dec。它不求解相机指向：切点必须由输入显式给出。整个
计算只用 Python 标准库，不访问网络、星表服务或天文图像。

```python
from star_crossmatch import solve_astrometry

result = solve_astrometry(data)
```

## 输入与拒绝边界

`data` 必须是对象，且字段恰为：

| 字段 | 规则 |
| --- | --- |
| `tangent_point` | 恰含 `ra`、`dec`，使用仓库十进制字符串语法；`ra∈[0,360)`、`dec∈[-90,90]` |
| `controls` | 1..1000 个对象；每行恰含 `id,x,y,ra,dec` |
| `detections` | 1..2000 个对象；每行恰含 `id,x,y` |
| `residual_limit_arcsec` | 最多 9 位整数、6 位小数的正十进制字符串 |
| `min_inliers` | JSON 整数（拒绝布尔值），范围 3..1000 |

ID 沿用 `[A-Za-z0-9_:-]{1,64}`；控制星 ID 必须在 `controls` 内唯一，检测 ID
必须在 `detections` 内唯一。两个身份命名空间相互独立，允许同一个字符串同时
作为控制星和检测 ID。`x,y` 沿用有界十进制字符串语法
`-?(?:0|[1-9][0-9]{0,8})(?:\.[0-9]{1,6})?`，不另设物理单位或视场上限。

缺字段、多余字段、错误类型、非有限值、数字/坐标越界、重复身份、空数组或
数组超界均抛 `ValueError`。切点或任一控制星不在切点所属的开可投影半球
（球心角 `<90°`；gnomonic 分母必须为正）也抛 `ValueError`。函数不修改传入
对象；所有验证先于几何结果返回，不会留下“半成功”的业务结果。

## 投影与仿射模型

角度计算使用 Python `math` 双精度浮点数。设切点为
`(α0,δ0)`，控制星为 `(α,δ)`，则：

```text
cos c = sin δ0 sin δ + cos δ0 cos δ cos(α-α0)
xi    = [cos δ sin(α-α0)] / cos c
eta   = [cos δ0 sin δ - sin δ0 cos δ cos(α-α0)] / cos c
```

输出中的 `xi`、`eta` 单位为度。仿射参数顺序固定为
`(a,b,c,d,e,f)`：

```text
xi  = a*x + b*y + c
eta = d*x + e*y + f
```

结果同时提供：

- `affine`：命名字典 `{a,b,c,d,e,f}`，也支持整数下标和解包为六元组；
- `affine_parameters`：按上述顺序排列的六个值；
- `matrix`：`[[a,b,c],[d,e,f]]`。

最小二乘在内部先对 `x,y` 中心化和尺度归一化，再解正规方程，最后换算回原
坐标单位。三个 F raction 精确共线检查用于秩判定，因此十进制控制点共线不会
被浮点误差误判为可解。

## 确定性稳健共识

1. 控制星先按 ID 排序；所有后续候选、评分和平局顺序均与输入行序无关。
2. 当三元组总数不超过 4096（即至多 30 个控制星）时，穷举所有 Fraction 精确
   非共线三元组。更大的输入使用固定种子 `0x51A57206` 的确定性
   `random.Random` 样本（至多 900 个不重复三元组）加最多 120 个按 ID 顺序的
   局部三元组；不使用系统随机源、时间或网络。
3. 每个候选三元组先拟合仿射解，再计算所有控制星的**球面角残差**。实现把
   切平面预测恢复为单位方向向量，用 `atan2(||u×v||,u·v)` 求角，残差单位为
   角秒；这不是只比较切平面线性残差。
4. `residual_arcsec <= residual_limit_arcsec` 的点为当前候选内点。随后用全部
   当前内点做最小二乘重拟合、重新认定内点，直到内点集合稳定（最多 8 轮）。
5. 只保留最终至少有 3 个内点的不同内点集合。合格解必须满足
   `inlier_count >= min_inliers`。
6. 合格解先按内点数降序，再按内点球面残差平方和
   （`residual_sum_squares_arcsec2`，单位 `arcsec²`）升序选择。SSE 仅允许
   `max(1e-20, 1e-12*|best_sse|)` 的双精度末位容差；该容差只用于识别数学上
   相同、由往返投影产生末位差异的平局，不改变残差上限本身的无 epsilon 比较。
7. 同内点数且 SSE 平局容差内的解均保留为最佳候选，并按排序后的内点 ID 元组
   排序。若任意两个最佳解对**任一输入检测**（不要求该检测在某个凸包内）给出
   的天球方向球面角差 `> residual_limit_arcsec`，整场状态为
   `ambiguous_solution`，不按行序或 ID 选解。边界等于上限不算超过。
8. 若所有最佳解对所有检测的位置差都不超过上限，则它们在给定科学容差下等价；
   此时选择内点 ID 元组字典序最小者，并在 `selection_rule` 中记录上述规则。

## 校准状态

顶层 `status` 与 `calibration.status` 的规范值为：

- `solved` / `calibrated`：找到唯一或在科学容差内等价的有效解；
- `rank_deficient`：全部控制星找不到三个非共线点（包括少于 3 个或所有有效
  点共线/重合）。此时不提供仿射解；
- `insufficient_inliers`：几何秩足够且搜索到了候选共识，但没有候选达到
  `min_inliers`；`inlier_count` 报告搜索到的最大共识规模；
- `ambiguous_solution`：存在两个或更多同分最佳解，并且它们会使至少一个检测
  的天球位置相差超过给定上限。`alternatives` 完整列出内点、SSE 和仿射参数。

顶层另有布尔 `valid`，仅 `solved` 时为 `true`。`insufficient_inliers` 时顶层
`affine` 仍为 `null`，但 `calibration.best_candidate` 按“最多内点、再最小 SSE、
再按 ID”披露搜索到的最佳不合格候选；对应控制星行给出该候选的残差和布尔内点
标记。为兼容常见命名，状态字符串与 `"ok"`、`"calibrated"`、`"ambiguous"`、
`"too_few_inliers"` 等近义写法的直接相等比较返回真，但 JSON 中的规范值仍是
上列字符串。

## 凸包与检测状态

有效解只对其**最终内点**的探测器坐标凸包负责。凸包用输入十进制对应的
`Fraction` 精确构造，边界包含在责任区内：

- 点在凸包内部、顶点或边上：用最终仿射参数计算 `xi,eta`，再做标准 gnomonic
  反投影；RA 用 `degrees(ra) % 360` 归一到 `[0,360)`；
- 点在凸包外：`status="outside_footprint"`，`ra`、`dec` 和 `coordinates`
  均为 `null`，不进行外推；
- 无有效解（秩不足、内点不足、歧义）：所有检测坐标为 `null`，检测状态使用
  对应的整场状态代码，`inside_footprint=null`。

`solved` 时内部/边界检测的规范状态是 `transformed`，并额外提供
`footprint_status="inside_footprint"` 与布尔 `inside_footprint=true`。

## 输出结构

```text
{
  "status": "solved",
  "calibration_status": "calibrated",
  "valid": true,
  "tangent_point": {"ra": "0", "dec": "0"},
  "residual_limit_arcsec": "0.05",
  "min_inliers": 6,
  "control_count": 10,
  "detection_count": 3,
  "inlier_count": 9,
  "affine": {"a": "...", "b": "...", "c": "...",
             "d": "...", "e": "...", "f": "..."},
  "affine_parameters": ["...", "...", "...", "...", "...", "..."],
  "matrix": [["...", "...", "..."], ["...", "...", "..."]],
  "residual_sum_squares_arcsec2": "...",
  "convex_hull": [["x", "y"], ...],
  "selection_rule": "固定、可复现的候选与排序说明",
  "calibration": {"status": "calibrated", ...},
  "controls": [
    {
      "id": "...", "x": "原始字符串", "y": "原始字符串",
      "ra": "原始字符串", "dec": "原始字符串",
      "xi": "浮点最短往返字符串（度）",
      "eta": "浮点最短往返字符串（度）",
      "residual_arcsec": "浮点最短往返字符串或 null",
      "residual": "同 residual_arcsec",
      "inlier": true,
      "is_inlier": true
    }
  ],
  "detections": [
    {
      "id": "...", "x": "原始字符串", "y": "原始字符串",
      "inside_footprint": true,
      "status": "transformed",
      "footprint_status": "inside_footprint",
      "ra": "浮点最短往返字符串",
      "dec": "浮点最短往返字符串",
      "coordinates": {"ra": "...", "dec": "..."}
    }
  ]
}
```

`controls` 与 `detections` 均按 ID Unicode/ASCII 字典序排序。秩不足时无法计算
仿射残差，因此 `affine` 为 `null`、控制星残差与内点状态为 `null`；内点不足时
控制星行反映 `calibration.best_candidate`；歧义时不存在唯一最终解，顶层
`residual_arcsec/inlier` 为 `null`，但每行的 `candidate_residuals` 与
`candidate_inliers` 按 `alternatives` 顺序披露每个并列解。

## 数值序列化

- 输入十进制字段保持其原始字符串并原样回显，不重新格式化。
- 所有三角/最小二乘计算均为 Python 双精度浮点数。
- `xi`、`eta`、输出 `ra/dec`、残差、SSE 和仿射系数用 `repr(float)` 生成最短
  往返十进制字符串；它们是 JSON 字符串，不声称精确有理值。该字符串子类仍可
  直接与数值做 `assertAlmostEqual`/算术比较。
- `id`、布尔值、整数计数和 `min_inliers` 保持 JSON 布尔/整数类型。
- 输出不含 `NaN` 或 `Infinity`；可直接通过 `star_crossmatch.io.publish` 写入
  JSON。`publish` 的原子新建、拒绝覆盖和 8 MiB 输入读取契约不变。

## 离线示例与测试

```sh
python3 -I examples/run_astrometry.py
python3 -I run_tests.py
```

示例生成原创合成仿射控制星、一个坏控制星和凸包内/边/外检测，写入
`run/astrometry.json`。`tests/test_astrometry.py` 覆盖：

1. 已知仿射解、小噪声和错误控制星；
2. 删除决定凸包边界的控制星后，检测由覆盖变为 `outside_footprint`；
3. RA 跨 0 度的 gnomonic 投影与反投影；
4. 共线控制点导致 `rank_deficient`；
5. 有效内点不足；
6. 两个等规模、同 SSE 且检测位置不同的最佳共识导致 `ambiguous_solution`；
7. 重排控制星/检测行后，JSON 结果按身份比较完全相同；
8. 字段、类型、十进制语法、坐标范围、身份重复、数组边界和半球拒绝。

这些测试不依赖网络或天文图像。
