# 整批视差贝叶斯距离推断契约（infer_distance_catalog）

新增 `infer_distance_catalog(data)`：把整批源的带符号视差与一个**全批共享**的加性视差零点 `z` 联合建模，返回每个源的距离后验分位数、零点后验摘要与指定源对的排序概率。仅依赖 Python 标准库，确定性、可复现，不访问网络。

## 概率模型

单位：视差与零点为 mas，距离为 pc，尺度长度 `scale_length_pc` 为 pc。

```
z      ~ Normal(mean_mas, sigma_mas^2)            # sigma_mas=0 退化为固定点质量
r_i    ~ p(r|L) = r^2 exp(-r/L) / (2 L^3),  r > 0
y_i | r_i, z ~ Normal(1000/r_i + z, sigma_i^2)    # 给定 (r,z) 条件独立
```

目标联合后验（忽略与 r、z 无关的归一化常数）：

```
p(z, r_1, …, r_n | y) ∝ p(z) · Π_i p(r_i|L) · Normal(y_i; 1000/r_i + z, sigma_i^2)
```

所有源通过同一个 `z` 耦合：逐源距离是该联合后验的**边缘分布**，零点摘要是 `z` 的全批边缘后验矩，pair 概率是对联合后验的二维边缘积分。任何源都不被当作固定距离/固定定标标准；`sigma_mas=0` 按固定零点先验的极限（点质量）处理。视差允许为零或负数，距离始终严格为正，全流程不使用 `1000/parallax` 作为估计或截断边界。

## 输入

顶层对象恰含四个字段：

- `sources`：1..64 行，每行恰含 `id`、`parallax_mas`、`sigma_mas`；`id` 沿用本仓 ASCII 身份语法且全批唯一；视差为可带符号的有界十进制字符串，`sigma_mas` 为严格正的十进制字符串。
- `zero_point`：恰含 `mean_mas`（可负）与 `sigma_mas`（非负，允许 `0`）的对象，值为有界十进制字符串。
- `scale_length_pc`：严格正的有界十进制字符串。
- `pairs`：0..256 行，每行恰含 `id`、`nearer`、`farther`；`id` 全批唯一，两个端点必须是已声明的不同 source id。字段名只声明身份，真实远近完全由后验概率决定。

布尔冒充数字、非字符串数值、非法十进制语法（如 `.5`、`1e3`、`nan`）、非有限值、空数组、超规模数组、重复身份、未知端点、相同端点、多余/缺失字段一律抛 `ValueError`；失败时不修改输入对象，也不返回部分结果。

## 数值方法（确定性）

距离维两次换元 `u = 1000/r`、`s = ln u`。先验 `g(u) ∝ u^-4 exp(-1000/(uL))`、高斯似然与雅可比合并后，`s` 空间条件后验为光滑单峰函数：

```
h(s) ∝ exp(-3s - 1000/(L e^s) - (y - z - e^s)^2 / (2 sigma^2))
```

- 对每个 `(源, z)`：确定性上坡 + 黄金分割定位众数，把积分区间扩张到相对对数密度低于 -46 的尾部，再用**全局误差预算的优先队列自适应 Simpson**（始终分裂误差最大的面板，相对容差 1e-11，面板数/最小宽度硬保护）求归一化常数；分位数由混合 CDF 在 `s` 空间固定步数二分反演后经 `r = 1000/e^s` 换回。
- 零点维：先做数据驱动的众数定位与二阶曲率估计，得到贴合后验形状的积分区间，再在其上用嵌套复合 Simpson 网格（32/64/128… 直到 4096 个等距区间，偶数节点位级复用）。权重归一化后给出零点后验矩；逐源边缘是各 `z` 条件后验按零点后验权重的混合。
- 排序概率：`r_n < r_f` 等价于 `u_n > u_f`（即 `s_n > s_f`），故 `P = Σ_z w_z ∫ p_n(s|z) (1 - F_f(s|z)) ds`，在近端 Simpson 面板上求积并按各条件自身质量归一。
- 所有相对密度都先减去对数峰值再取指数，避免直接指数下溢；非有限中间量、非正归一化常数或达到积分上限仍不收敛都抛 `ValueError`。

距离分位数在相邻两级零点网格的最大变化低于 8e-5 pc 时接受（严格优于交付的 1e-4 pc）；pair 概率的误差由距离维 Simpson 单独保证在 1e-5 以内。持久化数值用紧凑 `array('d')`，并按当前零点网格原地淘汰上一级条件后验以控制内存。

该引擎（`s=ln(1000/r)` 空间自适应 Simpson + 零点复合 Simpson）与测试侧独立 oracle（`s=ln r` 空间自适应 Gauss-Legendre + 零点梯形网格，见 `tests/_distance_oracle.py`）是两套独立实现，互不导入、不共用积分代码。

## 输出

可直接 `json.dumps(..., allow_nan=False)` 的版本化 dict，`format` 固定为 `star-crossmatch-distance-catalog/1`：

- `zero_point`：`{mean_mas, sigma_mas}`，共享零点的后验均值与标准差（固定零点先验时 `sigma_mas` 为 `0.0`）。
- `sources`：按 source id 字典序，每行 `{id, median_pc, ci_lower_pc, ci_upper_pc}`；下界/中位/上界对应等尾 68.27% 区间（分位点 0.15865、0.5、0.84135），均为有限、严格为正且严格有序的 pc 数值。
- `pairs`：按 pair id 字典序，每行 `{id, nearer, farther, probability}`，`probability` 为 `[0,1]` 内的有限数；交换 nearer/farther 的两个互补 pair 概率之和在数值容差内为 1。

重排输入的 source 行或 pair 行不改变按身份比较的 JSON 结果。固定小批次的独立高精度参考值保存在 `tests/distance_oracle_fixture.json`，离线测试以距离 1e-4 pc、pair 1e-5 的绝对容差对照。
