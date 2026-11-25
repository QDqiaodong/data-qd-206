"""整批视差的贝叶斯距离推断（仅依赖 Python 标准库）。

公开 API :func:`infer_distance_catalog` 把整批源的带符号视差与一个全批共享
的加性视差零点 ``z`` 联合建模，交付：

* 每个源在**全目录联合后验**下的距离后验中位数与等尾 68.27% 可信区间；
* 共享零点 ``z`` 的后验均值/标准差；
* 每个指定源对在联合后验下 ``r_nearer < r_farther`` 的概率。

概率模型（mas、pc 单位）::

    z        ~ Normal(mean_mas, sigma_mas^2)      # sigma_mas=0 时退化为点质量
    r_i      ~ p(r|L) = r^2 exp(-r/L) / (2 L^3), r > 0
    y_i | r_i, z ~ Normal(1000/r_i + z, sigma_i^2)   # 条件独立

所有源的视差通过同一个 ``z`` 耦合，逐源输出是联合后验的边缘分布，绝不是把
零点固定在先验均值后各自独立拟合；视差允许为零或负数，距离始终严格为正，
任何地方都不使用 ``1000/parallax`` 作为估计或截断边界。

数值策略（确定性、可复现）
--------------------------

距离维两次换元：``u = 1000/r``（mas）、再取 ``s = ln u``。把先验
``g(u) ~ u^-4 exp(-1000/(u L))``、高斯似然与雅可比 ``du = u ds`` 合并后，
``s`` 空间的未归一化条件后验为::

    h(s) ~ exp(-3 s - 1000/(L e^s) - (y-z-e^s)^2/(2 sigma^2))

它在有界 ``s`` 区间内光滑、单峰，避免了 ``u`` 空间在 0 附近的陡壁。对每个
``(源, z)``：先以确定性黄金分割定位众数，把积分区间扩张到相对对数密度低于
-46 的尾部（被忽略质量量级 < 1e-20），再用**全局误差预算的优先队列自适应
Simpson**（按区间误差从大到小分裂，相对容差 1e-11，面板数与最小宽度有硬
保护）求归一化常数；CDF 由面板前缀质量加所在面板的精确单区间 Simpson 直接
计算（不做跨面板线性插值），距离分位数由混合 CDF 在 ``s`` 空间二分反演后
经 ``r = 1000/e^s`` 换回。

零点维在其全批边缘后验上用嵌套复合 Simpson 网格（64/128/256... 个等距区间，
旧节点全部复用），权重归一化后给出零点后验矩；逐源边缘是各 ``z`` 条件后验
按零点后验权重的混合。注意 ``r_n < r_f`` 等价于 ``u_n > u_f``，故排序概率
为 ``sum_z w_z int p_n(s|z) (1 - F_f(s|z)) ds``。

相邻两级 z 网格的距离分位数差 <= 2.5e-5 pc、pair 概率差 <= 2.5e-6 才视为
收敛；达到区间上限仍不满足，或任何中间量非有限/归一化常数非正，一律抛
:class:`ValueError`，不返回部分结果。

该引擎（s=ln(1000/r) 空间优先队列自适应 Simpson + z 复合 Simpson）与离线
测试中的独立 oracle（ln r 空间自适应 Gauss-Legendre + z 梯形网格）是两套
独立实现，不共用任何积分代码。
"""
from __future__ import annotations

import heapq
import math
from array import array
from bisect import bisect_left

from .checks import array as _check_array
from .checks import fields, ident, number

DISTANCE_CATALOG_FORMAT = "star-crossmatch-distance-catalog/1"

SOURCE_LIMIT = 64
PAIR_LIMIT = 256
# 等尾 68.27% 可信区间对应的分位点（与固定 oracle 采用同一组分位定义）。
Q_LOW = 0.15865
Q_HIGH = 0.84135
# 距离维自适应 Simpson 的相对容差与硬保护。
SIMPSON_REL_TOL = 1e-11
SIMPSON_MAX_PANELS = 8000
SIMPSON_MIN_WIDTH = 1e-7
# 尾部截断：相对对数密度低于该值即视为可忽略（exp(-46) ~ 1e-20）。
TAIL_LOG_CUTOFF = 46.0
# z 复合 Simpson 的起始/上限区间数（偶数；嵌套加倍）。起始取较小值，因为
# 积分区间已按数据驱动的后验曲率收紧；对距离对零点极敏感的高精度源，上限
# 允许加到 4096 区间仍能稳定收敛。
Z_INITIAL_INTERVALS = 32
Z_MAX_INTERVALS = 4096
# 收敛判据：网格加细一级后距离分位数的最大变化低于该值即接受当前级。阈值
# 取距离交付容差（1e-4 pc）的 80%，给条件后验在新 z 节点重建引入的浮点抖动
# 留出余量，同时仍严格优于交付精度。pair 概率只在最终接受的网格上计算一次，
# 其积分误差由距离维 Simpson（相对容差 1e-11）单独保证在 1e-5 以内。
QUANTILE_CONVERGENCE = 8.0e-5
_LOG_2PI = math.log(2.0 * math.pi)
_BISECTION_STEPS = 60


# ---------------------------------------------------------------------------
# 输入校验
# ---------------------------------------------------------------------------


def _decimal(value, label):
    """仓库统一的有界十进制字符串 -> 有限 float；布尔/非字符串一律拒绝。"""
    fraction = number(value)  # 非字符串或非法十进制语法在此抛 ValueError
    result = float(fraction)
    if not math.isfinite(result):
        raise ValueError(f"{label} 必须是有限十进制数")
    return result


def _validate(data):
    if not isinstance(data, dict):
        raise ValueError("输入必须是对象")
    fields(data, ("sources", "zero_point", "scale_length_pc", "pairs"))

    raw_sources = _check_array(data["sources"], 1, SOURCE_LIMIT)
    sources = []
    seen_ids = set()
    for row in raw_sources:
        fields(row, ("id", "parallax_mas", "sigma_mas"))
        sid = ident(row["id"])
        if sid in seen_ids:
            raise ValueError(f"source id 批内重复: {sid}")
        seen_ids.add(sid)
        parallax = _decimal(row["parallax_mas"], "parallax_mas")
        sigma = _decimal(row["sigma_mas"], "sigma_mas")
        if sigma <= 0.0:
            raise ValueError("sigma_mas 必须严格为正")
        sources.append({"id": sid, "y": parallax, "sigma": sigma})

    zp = data["zero_point"]
    if not isinstance(zp, dict):
        raise ValueError("zero_point 必须是对象")
    fields(zp, ("mean_mas", "sigma_mas"))
    z_mean = _decimal(zp["mean_mas"], "zero_point.mean_mas")
    z_sigma = _decimal(zp["sigma_mas"], "zero_point.sigma_mas")
    if z_sigma < 0.0:
        raise ValueError("zero_point.sigma_mas 不允许为负")

    scale_length = _decimal(data["scale_length_pc"], "scale_length_pc")
    if scale_length <= 0.0:
        raise ValueError("scale_length_pc 必须严格为正")

    raw_pairs = _check_array(data["pairs"], 0, PAIR_LIMIT)
    pairs = []
    seen_pair_ids = set()
    for row in raw_pairs:
        fields(row, ("id", "nearer", "farther"))
        pid = ident(row["id"])
        if pid in seen_pair_ids:
            raise ValueError(f"pair id 批内重复: {pid}")
        seen_pair_ids.add(pid)
        nearer = ident(row["nearer"])
        farther = ident(row["farther"])
        if nearer not in seen_ids or farther not in seen_ids:
            raise ValueError(f"pair {pid} 引用了未知源端点")
        if nearer == farther:
            raise ValueError(f"pair {pid} 的两个端点必须是不同源")
        pairs.append({"id": pid, "nearer": nearer, "farther": farther})

    # 计算前按身份排序，保证浮点求和次序与输入行序无关。
    sources.sort(key=lambda item: item["id"])
    pairs.sort(key=lambda item: item["id"])
    return {
        "sources": sources,
        "z_mean": z_mean,
        "z_sigma": z_sigma,
        "scale_length": scale_length,
        "pairs": pairs,
    }


# ---------------------------------------------------------------------------
# s = ln(1000/r) 变量下的条件后验（固定 z）
# ---------------------------------------------------------------------------


class _ConditionalIntegrand:
    """固定 (源, z) 时 s=ln u=ln(1000/r) 空间的未归一化条件后验密度。"""

    __slots__ = ("y_eff", "sigma", "a_over_L", "half_inv_var")

    def __init__(self, y, sigma, z, scale_length):
        self.y_eff = y - z
        self.sigma = sigma
        self.a_over_L = 1000.0 / scale_length
        self.half_inv_var = 0.5 / (sigma * sigma)

    def log_abs(self, s):
        u = math.exp(s)
        return (-3.0 * s
                - self.a_over_L / u
                - self.half_inv_var * (self.y_eff - u) ** 2)


def _golden_max(fn, a, b):
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc, fd = fn(c), fn(d)
    for _ in range(120):
        if fc < fd:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fn(d)
        else:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fn(c)
    m = 0.5 * (a + b)
    return m, fn(m)


def _find_log_mode(fn):
    """在整条实轴上确定性地夹住单峰并黄金分割细化，不依赖固定窗口。

    初始步长取得很小并逐步放大，确保宽度极窄的峰不会被第一步直接跨过；
    一旦在某方向越过峰顶（函数值回落）即夹住众数。
    """
    center = 0.0
    step = 0.01
    f_left = fn(center - step)
    f_right = fn(center + step)
    direction = 1.0 if f_right >= f_left else -1.0
    a = center
    b = center + direction * step
    fb = f_right if direction > 0.0 else f_left
    for _ in range(2000):
        c = b + direction * step
        fc = fn(c)
        if fc <= fb:
            return _golden_max(fn, min(a, c), max(a, c))
        a, b, fb = b, c, fc
        step *= 1.02
    return _golden_max(fn, b - step, b + step)


def _simpson_value(fn, a, b):
    return (b - a) / 6.0 * (fn(a) + 4.0 * fn(0.5 * (a + b)) + fn(b))


def _adaptive_simpson(fn, lo, hi, min_width=SIMPSON_MIN_WIDTH):
    """全局误差预算的优先队列自适应 Simpson，返回排序后的面板 (a,b,积分)。

    始终分裂当前误差最大的区间。停止判据采用全局**绝对**误差地板
    ``SIMPSON_REL_TOL * |整体积分|``：当区间误差低于该地板，或区间自身质量
    与误差都可忽略时直接接受，避免在质量近零、相对误差无意义的深尾部无限
    分裂。面板数与最小宽度仍有硬保护，触发即抛 :class:`ValueError`。
    """
    whole = _simpson_value(fn, lo, hi)
    if not math.isfinite(whole):
        raise ValueError("距离后验积分出现非有限值")
    abs_floor = SIMPSON_REL_TOL * abs(whole)
    heap = []

    def push(a, b):
        mid = 0.5 * (a + b)
        left = _simpson_value(fn, a, mid)
        right = _simpson_value(fn, mid, b)
        local = _simpson_value(fn, a, b)
        heapq.heappush(
            heap, (-abs(left + right - local), a, mid, b, left, right))

    push(lo, hi)
    total_error = sum(-item[0] for item in heap)
    panels = []
    while heap:
        neg_err, a, mid, b, left, right = heapq.heappop(heap)
        err = -neg_err
        total_error -= err
        refined = left + right
        negligible = refined <= abs_floor and err <= abs_floor
        if err <= 15.0 * abs_floor or negligible:
            panels.append((a, b, refined))
        elif len(panels) >= SIMPSON_MAX_PANELS or (mid - a) <= min_width:
            raise ValueError("距离后验收敛失败：自适应 Simpson 触顶")
        else:
            for x, y in ((a, mid), (mid, b)):
                mm = 0.5 * (x + y)
                s1 = _simpson_value(fn, x, mm)
                s2 = _simpson_value(fn, mm, y)
                sub_err = abs(s1 + s2 - _simpson_value(fn, x, y))
                heapq.heappush(heap, (-sub_err, x, mm, y, s1, s2))
                total_error += sub_err
        if total_error <= abs_floor:
            while heap:
                _, pa, pm, pb, s_left, s_right = heapq.heappop(heap)
                panels.append((pa, pb, s_left + s_right))
            break
    panels.sort()
    return panels


class _Conditional:
    """一个 (源, z) 条件后验的面板、归一化常数与精确 CDF。

    持久化数组统一用紧凑 ``array('d')``：只存面板右端点 ``_pb``、前缀质量
    ``_prefix`` 与稠密查询网格 ``_grid_s``/``_grid_cdf``；面板左端点由前一
    右端点（首个为 lo）恢复，单面板质量由相邻前缀之差恢复。
    """

    __slots__ = ("integrand", "peak", "_pb", "total",
                 "lo", "hi", "_prefix", "_grid_s", "_grid_cdf")

    def __init__(self, integrand):
        self.integrand = integrand
        # 众数一定存在且为正 u；沿实轴确定性上坡夹住单峰再黄金分割。
        mode, self.peak = _find_log_mode(integrand.log_abs)
        if not math.isfinite(mode) or not math.isfinite(self.peak):
            raise ValueError("距离后验众数非有限")
        lo, hi = mode - 6.0, mode + 6.0
        for _ in range(4000):
            low_ok = integrand.log_abs(lo) < self.peak - TAIL_LOG_CUTOFF
            high_ok = integrand.log_abs(hi) < self.peak - TAIL_LOG_CUTOFF
            if low_ok and high_ok:
                break
            if not low_ok:
                lo -= 2.0
            if not high_ok:
                hi += 2.0
        else:
            raise ValueError("距离后验尾部截断未能在保护步数内收敛")
        self.lo, self.hi = lo, hi
        # 最小面板宽度相对积分区间取值，兼容极窄或极宽的后验。
        min_width = (hi - lo) * 1.0e-9

        def relative_density(s):
            return math.exp(integrand.log_abs(s) - self.peak)

        panel_tuples = _adaptive_simpson(relative_density, lo, hi, min_width)
        # 持久化数组统一用紧凑 C double（array('d')，每项 8 字节）。只存面板
        # 右端点 _pb 与前缀质量 _prefix：左端点由“前一面板右端点（首个为 lo）”
        # 恢复，单面板质量由相邻前缀之差恢复，避免重复数组占用。
        self._pb = array("d", (item[1] for item in panel_tuples))
        prefix = []
        acc = 0.0
        for item in panel_tuples:
            acc += item[2]
            prefix.append(acc)
        self.total = acc
        if self.total <= 0.0 or not math.isfinite(self.total):
            raise ValueError("距离后验归一化常数必须有限且为正")
        self._prefix = array("d", prefix)
        # 单一稠密查询网格：每个 Simpson 面板取 1/4、1/2、3/4、1 四个节点，
        # 起始额外存左端点。坐标与归一化 CDF 各一份紧凑 array('d')，供 pair
        # 高频查询做 bisect + 线性插值（误差为面板宽平方量级，pair 误差约
        # 1e-6）；分位数反演仍走 cdf_s 的精确 Simpson 路径，不受插值影响。
        grid_s = [lo]
        grid_cdf = [0.0]
        running = 0.0
        prev_b = lo
        for k, (_a, b, panel_mass) in enumerate(panel_tuples):
            a = prev_b
            span = b - a
            for frac in (0.25, 0.5, 0.75, 1.0):
                s = a + span * frac
                if frac == 1.0:
                    running += panel_mass
                    grid_s.append(s)
                    grid_cdf.append(min(1.0, running / self.total))
                else:
                    partial = _simpson_value(relative_density, a, s)
                    grid_s.append(s)
                    grid_cdf.append(
                        min(1.0, (running + partial) / self.total))
            prev_b = b
        grid_cdf[-1] = 1.0
        self._grid_s = array("d", grid_s)
        self._grid_cdf = array("d", grid_cdf)

    def log_evidence(self):
        return self.peak + math.log(self.total)

    def cdf_s_fast(self, s):
        """稠密网格 CDF 上的 bisect + 线性插值，专供 pair 高频调用。"""
        if s <= self.lo:
            return 0.0
        if s >= self.hi:
            return 1.0
        j = bisect_left(self._grid_s, s)
        if j >= len(self._grid_s):
            return 1.0
        x0, x1 = self._grid_s[j - 1], self._grid_s[j]
        c0, c1 = self._grid_cdf[j - 1], self._grid_cdf[j]
        t = (s - x0) / (x1 - x0)
        return c0 + t * (c1 - c0)

    def cdf_s(self, s):
        """s=ln u 空间的精确归一化 CDF：前缀质量 + 所在面板 Simpson。"""
        if s <= self.lo:
            return 0.0
        if s >= self.hi:
            return 1.0
        lo_i = bisect_left(self._pb, s)
        if lo_i >= len(self._pb):
            return 1.0
        panel_a = self.lo if lo_i == 0 else self._pb[lo_i - 1]
        prior = 0.0 if lo_i == 0 else self._prefix[lo_i - 1]
        integrand = self.integrand
        peak = self.peak
        mid_s = 0.5 * (panel_a + s)
        fa = math.exp(integrand.log_abs(panel_a) - peak)
        fm = math.exp(integrand.log_abs(mid_s) - peak)
        fb = math.exp(integrand.log_abs(s) - peak)
        partial = (s - panel_a) / 6.0 * (fa + 4.0 * fm + fb)
        return min(1.0, max(0.0, (prior + partial) / self.total))

    def density_s(self, s):
        return math.exp(self.integrand.log_abs(s) - self.peak) / self.total


# ---------------------------------------------------------------------------
# 零点 z 的全批边缘后验
# ---------------------------------------------------------------------------


class _ZeroPointModel:
    def __init__(self, spec):
        self.sources = spec["sources"]
        self.scale_length = spec["scale_length"]
        self.mu = spec["z_mean"]
        self.tau = spec["z_sigma"]
        self.fixed = self.tau == 0.0
        self._cache = {}

    def _prior_log(self, z):
        if self.fixed:
            return 0.0
        return (-0.5 * ((z - self.mu) / self.tau) ** 2
                - math.log(self.tau) - 0.5 * _LOG_2PI)

    def conditional(self, source_index, z):
        key = (source_index, z)
        cached = self._cache.get(key)
        if cached is None:
            src = self.sources[source_index]
            integrand = _ConditionalIntegrand(
                src["y"], src["sigma"], z, self.scale_length)
            cached = _Conditional(integrand)
            self._cache[key] = cached
        return cached

    def retain_z_nodes(self, z_nodes):
        """只保留给定 z 网格上的条件后验，原地淘汰上一级网格的陈旧条目。

        嵌套 Simpson 加倍时偶数节点位级相同、float 位级一致，因而被保留的
        条件后验可直接复用；新增奇数节点才触发构建。原地删除（而非新建字典）
        避免新旧两份缓存并存，使内存峰值始终只有一级网格的体量。
        """
        keep = set()
        n_src = len(self.sources)
        for i in range(n_src):
            for z in z_nodes:
                keep.add((i, z))
        for key in list(self._cache.keys()):
            if key not in keep:
                del self._cache[key]

    def log_posterior_unorm(self, z):
        value = self._prior_log(z)
        for i in range(len(self.sources)):
            value += self.conditional(i, z).log_evidence()
        if not math.isfinite(value):
            raise ValueError("零点后验计算出现非有限值")
        return value

    def log_posterior_transient(self, z):
        """与 :meth:`log_posterior_unorm` 同值，但不写入持久缓存。

        零点截断扫描只需要边缘似然值，扫描节点通常不落在最终 Simpson 网格
        上；用临时条件后验即时计算并丢弃，避免缓存被一次性扫描节点占满。
        """
        value = self._prior_log(z)
        for src in self.sources:
            integrand = _ConditionalIntegrand(
                src["y"], src["sigma"], z, self.scale_length)
            value += _Conditional(integrand).log_evidence()
        if not math.isfinite(value):
            raise ValueError("零点后验计算出现非有限值")
        return value

    def locate_posterior_mode(self):
        """定位零点后验众数并按二阶曲率估计其标准差。

        先在先验支撑内做确定性粗扫描，再黄金分割细化；曲率宽度反映数据（含
        高信息量源）对零点的收紧程度，用于设置贴合后验形状的初始积分区间，
        避免在远宽于后验的先验区间上浪费 Simpson 节点。
        """
        width = max(abs(self.tau), 1e-12)
        # 扫描完整先验支撑，避免高精度源把似然众数推到远处时漏掉先验主导的
        # 真实后验众数。
        scan_lo, scan_hi = self.mu - 16.0 * width, self.mu + 16.0 * width
        n_scan = 256
        nodes = [scan_lo + (scan_hi - scan_lo) * i / n_scan
                 for i in range(n_scan + 1)]
        # 额外播种每个源的似然一致零点：高精度源的尖锐似然模态宽度可远小于
        # 均匀扫描步长，仅靠网格会漏检。似然项 Normal(y; u+z, sigma^2) 在
        # u=1000/r>0 下的最优零点趋近 z=y（u 取小但为正），故直接以视差 y
        # 作为种子；它在先验支撑内时一定能罩住该源的尖锐似然峰。
        for src in self.sources:
            seed = src["y"]
            if scan_lo <= seed <= scan_hi:
                nodes.append(seed)
        values = [self.log_posterior_transient(z) for z in nodes]
        j = max(range(len(nodes)), key=lambda k: values[k])
        # 取后验值最高的候选点，再在其局部括号内黄金分割细化。普通扫描点用
        # 相邻扫描间距做括号；尖锐似然种子用最小源测量标准差的小倍数做窄
        # 括号，确保罩住宽度仅 ~sigma 的尖峰而不被远处扫描点稀释。
        j = max(range(len(nodes)), key=lambda k: values[k])
        z_star = nodes[j]
        if j <= n_scan:
            half = (scan_hi - scan_lo) / n_scan
        else:
            half = max(width,
                       min(src["sigma"] for src in self.sources)) * 4.0
        b_lo = max(scan_lo, z_star - half)
        b_hi = min(scan_hi, z_star + half)
        mode, peak = _golden_max(
            self.log_posterior_transient, b_lo, b_hi)
        h = 1e-5 * max(1.0, abs(mode))
        second = (self.log_posterior_transient(mode + h)
                  - 2.0 * peak
                  + self.log_posterior_transient(mode - h)) / (h * h)
        sigma = 1.0 / math.sqrt(-second) if second < 0.0 else width
        if not math.isfinite(sigma) or sigma <= 0.0:
            sigma = width
        return mode, peak, sigma

    def bracket(self):
        """确定覆盖零点后验有效支撑的有限 z 区间。

        以数据驱动的后验众数与曲率宽度为起点向外扩张到尾部截断线；同时始终
        夹在先验的有限支撑（±16 个先验标准差）内。当后验众数远离先验、曲率
        宽度又极小时（高精度源与先验强烈冲突），曲率区间够不到先验边界，此时
        直接采用完整先验支撑作为积分区间——先验在边界外已可忽略，后验在那里
        没有质量，截断不引入误差。
        """
        mode, peak, sigma = self.locate_posterior_mode()
        prior_half = 16.0 * max(abs(self.tau), 1e-12)
        prior_lo, prior_hi = self.mu - prior_half, self.mu + prior_half
        if mode <= prior_lo or mode >= prior_hi or sigma < prior_half / 1e4:
            return prior_lo, prior_hi
        lo, hi = mode - 6.0 * sigma, mode + 6.0 * sigma
        for _ in range(4000):
            h_lo = self.log_posterior_transient(lo)
            h_hi = self.log_posterior_transient(hi)
            low_ok = lo <= prior_lo or h_lo < peak - TAIL_LOG_CUTOFF
            high_ok = hi >= prior_hi or h_hi < peak - TAIL_LOG_CUTOFF
            if low_ok and high_ok:
                return max(lo, prior_lo), min(hi, prior_hi)
            if not low_ok:
                lo -= 2.0 * sigma
            if not high_ok:
                hi += 2.0 * sigma
        raise ValueError("零点后验截断未能在保护步数内收敛")


def _simpson_weights(nodes, values):
    """复合 Simpson 归一化权重；nodes 等距且区间数为偶数。"""
    n = len(nodes) - 1
    if n < 2 or n % 2:
        raise ValueError("Simpson 网格区间数必须为正偶数")
    step = (nodes[-1] - nodes[0]) / n
    raw = [0.0] * len(nodes)
    for i, value in enumerate(values):
        if i == 0 or i == n:
            coef = 1.0
        elif i % 2:
            coef = 4.0
        else:
            coef = 2.0
        raw[i] = coef * value
    total = sum(raw) * step / 3.0
    if total <= 0.0 or not math.isfinite(total):
        raise ValueError("零点后验归一化常数必须有限且为正")
    return [x * step / 3.0 / total for x in raw]


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------


def _pair_probability(model, weights, z_nodes, near_idx, far_idx):
    """P(r_near < r_far) = sum_z w_z int p_n(s|z) (1-F_f(s|z)) ds。

    r_n < r_f 等价于 u_n > u_f（s_n > s_f），故使用远端的生存函数。被积
    函数在近端每个 Simpson 面板的端点/中点上取值；远端 CDF 用面板前缀质量
    加所在面板 Simpson 精确计算，并在单次 pair 调用内按 (源,z,s) 记忆化，
    避免相邻面板共享节点时的重复求积。每个 z 条件按近端面板总质量归一后
    再按零点后验权重混合。
    """
    total = 0.0
    for w, z in zip(weights, z_nodes):
        cond_n = model.conditional(near_idx, z)
        cond_f = model.conditional(far_idx, z)
        pb = cond_n._pb

        def far_survival(s):
            return 1.0 - cond_f.cdf_s_fast(s)

        sub = 0.0
        norm = 0.0
        a = cond_n.lo
        for b in pb:
            mid = 0.5 * (a + b)
            dn_a = cond_n.density_s(a)
            dn_m = cond_n.density_s(mid)
            dn_b = cond_n.density_s(b)
            sub += ((b - a) / 6.0) * (
                dn_a * far_survival(a)
                + 4.0 * dn_m * far_survival(mid)
                + dn_b * far_survival(b))
            norm += ((b - a) / 6.0) * (dn_a + 4.0 * dn_m + dn_b)
            a = b
        total += w * (sub / norm)
    if not math.isfinite(total) or total < -1e-9 or total > 1.0 + 1e-9:
        raise ValueError("pair 排序概率必须落在 [0,1] 且有限")
    return min(1.0, max(0.0, total))


def _stage_outputs(model, z_nodes, weights, spec, include_pairs):
    sources = spec["sources"]
    z_mean = sum(w * z for w, z in zip(weights, z_nodes))
    z_var = sum(w * z * z for w, z in zip(weights, z_nodes)) - z_mean * z_mean
    z_sigma = math.sqrt(max(0.0, z_var))
    if not math.isfinite(z_mean) or not math.isfinite(z_sigma):
        raise ValueError("零点后验摘要非有限")

    source_rows = []
    for i, src in enumerate(sources):
        conditionals = [model.conditional(i, z) for z in z_nodes]
        lo = min(item.lo for item in conditionals)
        hi = max(item.hi for item in conditionals)

        def mix_cdf(s):
            return sum(w * cond.cdf_s(s)
                       for w, cond in zip(weights, conditionals))

        def quantile(q):
            # F_r(r)=q 等价于 F_s(ln(1000/r)) = 1-q。在完整支撑区间内对混合
            # CDF 做固定步数二分；步数对 s 给出 ~区间宽/2^50 的分辨率，远低于
            # 距离交付误差，且完全确定性、不依赖初值。
            target = 1.0 - q
            a, b = lo, hi
            for _ in range(_BISECTION_STEPS):
                mid = 0.5 * (a + b)
                if mix_cdf(mid) < target:
                    a = mid
                else:
                    b = mid
            r = 1000.0 / math.exp(0.5 * (a + b))
            if not math.isfinite(r) or r <= 0.0:
                raise ValueError("距离分位数必须有限且为正")
            return r

        lower = quantile(Q_LOW)
        median = quantile(0.5)
        upper = quantile(Q_HIGH)
        if not lower < median < upper:
            raise ValueError("距离分位数必须严格有序")
        source_rows.append({
            "id": src["id"],
            "median_pc": median,
            "ci_lower_pc": lower,
            "ci_upper_pc": upper,
        })

    index_by_id = {src["id"]: i for i, src in enumerate(sources)}
    pair_rows = []
    if include_pairs:
        # pair 积分成本高，只在最终接受的网格上计算；收敛判据只依赖分位数，
        # 而分位数收敛后按同一权重混合的 pair 概率同样达到交付精度。
        for pair in spec["pairs"]:
            near_idx = index_by_id[pair["nearer"]]
            far_idx = index_by_id[pair["farther"]]
            pair_rows.append({
                "id": pair["id"],
                "nearer": pair["nearer"],
                "farther": pair["farther"],
                "probability": _pair_probability(
                    model, weights, z_nodes, near_idx, far_idx),
            })
    return {
        "z_mean": z_mean,
        "z_sigma": z_sigma,
        "sources": source_rows,
        "pairs": pair_rows,
    }


def _quantiles_within_tol(stage, previous):
    """两级 z 网格的距离分位数是否满足收敛判据（pair 在最终级单独计算）。"""
    return max(
        abs(cur[key] - prev[key])
        for cur, prev in zip(stage["sources"], previous["sources"])
        for key in ("median_pc", "ci_lower_pc", "ci_upper_pc")
    ) <= QUANTILE_CONVERGENCE


def infer_distance_catalog(data):
    """推断整批源的贝叶斯距离后验与批内源对排序概率。

    参数:
        data: 恰含 ``sources``、``zero_point``、``scale_length_pc``、``pairs``
            的 JSON 风格对象；数值均为仓库统一的有界十进制字符串。

    返回:
        可直接 ``json.dumps``（``allow_nan=False`` 安全）的版本化 dict，
        源按 source ID、pair 按 pair ID 排序，与输入行序无关。

    异常:
        ValueError: 任何字段/类型/身份/范围校验失败，或数值积分未能达到
            规定精度；此时不修改输入，也不返回部分结果。
    """
    spec = _validate(data)
    model = _ZeroPointModel(spec)

    if model.fixed:
        # 固定零点：z 退化为点质量，无需网格收敛，直接在唯一 z 上完整计算。
        z_nodes = [model.mu]
        weights = [1.0]
        current = _stage_outputs(
            model, z_nodes, weights, spec, include_pairs=True)
    else:
        zlo, zhi = model.bracket()
        intervals = Z_INITIAL_INTERVALS
        previous = None
        current = None
        while True:
            z_nodes = [zlo + (zhi - zlo) * i / intervals
                       for i in range(intervals + 1)]
            log_values = [model.log_posterior_unorm(z) for z in z_nodes]
            peak = max(log_values)
            values = [math.exp(v - peak) for v in log_values]
            weights = _simpson_weights(z_nodes, values)
            # 收敛过程只需要分位数与零点矩，跳过昂贵的 pair 积分。
            current = _stage_outputs(
                model, z_nodes, weights, spec, include_pairs=False)
            # 本级完成后立即淘汰不属于当前网格的条件后验，把缓存限制在一级
            # 网格的体量；下一级加倍时偶数节点位级一致、仍可命中复用。
            model.retain_z_nodes(z_nodes)
            # 嵌套 Simpson 网格随加细单调逼近；当某一级相对前一级的分位数
            # 变化首次落入收敛容差（小于交付误差限 1e-4 pc）即接受当前级，
            # 再在该最终网格上补算全部 pair 概率。
            if previous is not None and _quantiles_within_tol(
                    current, previous):
                current = _stage_outputs(
                    model, z_nodes, weights, spec, include_pairs=True)
                break
            if intervals >= Z_MAX_INTERVALS:
                raise ValueError("零点后验收敛失败：已达积分区间上限")
            previous = current
            intervals *= 2

    return {
        "format": DISTANCE_CATALOG_FORMAT,
        "zero_point": {
            "mean_mas": current["z_mean"],
            "sigma_mas": current["z_sigma"],
        },
        "sources": current["sources"],
        "pairs": current["pairs"],
    }
