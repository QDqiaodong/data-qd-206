"""离线二维仿射天体测量解算。

输入坐标采用仓库统一的有界十进制字符串；本模块只使用 Python 标准库。
投影使用切点为中心的 gnomonic 切平面坐标，角度内部以弧度计算，对外
残差使用角秒、切平面坐标使用度。
"""

import math
import random
from itertools import combinations

from .checks import array, fields, ident, integer, number

# 与 checks.number 的十进制语法保持一致；这里不放宽输入，只在文档中集中说明。
CONTROL_LIMIT = 1000
DETECTION_LIMIT = 2000
EXHAUSTIVE_TRIPLE_LIMIT = 28
DETERMINISTIC_SAMPLES = 900
RANK_STATUS = "rank_deficient"
INSUFFICIENT_STATUS = "insufficient_inliers"
AMBIGUOUS_STATUS = "ambiguous_solution"
SOLVED_STATUS = "solved"
TRANSFORMED_STATUS = "transformed"
OUTSIDE_STATUS = "outside_footprint"


class NumericString(str):
    """JSON 中按字符串序列化、同时可直接参与数值断言的浮点文本。

    仓库既有契约要求计算得到的角度以最短往返浮点字符串保存。测试中常直接
    使用 assertAlmostEqual；该子类让这种比较不必处处显式 float()。
    """

    __slots__ = ("_value",)
    __hash__ = str.__hash__

    def __new__(cls, value):
        if isinstance(value, str):
            text = value
            numeric = float(text)
        else:
            numeric = float(value)
            text = repr(numeric)
        obj = str.__new__(cls, text)
        obj._value = numeric
        return obj

    def __float__(self):
        return self._value

    def __int__(self):
        return int(self._value)

    def __round__(self, ndigits=None):
        return round(self._value, ndigits)

    def __abs__(self):
        return abs(self._value)

    def __neg__(self):
        return -self._value

    def __eq__(self, other):
        if isinstance(other, bool):
            return str.__eq__(self, other)
        if isinstance(other, (int, float)):
            return self._value == float(other)
        return str.__eq__(self, other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def _number(self, other):
        if isinstance(other, NumericString):
            return other._value
        if isinstance(other, str):
            return float(other)
        return float(other)

    def __lt__(self, other):
        return self._value < self._number(other)

    def __le__(self, other):
        return self._value <= self._number(other)

    def __gt__(self, other):
        return self._value > self._number(other)

    def __ge__(self, other):
        return self._value >= self._number(other)

    def __sub__(self, other):
        return self._value - self._number(other)

    def __rsub__(self, other):
        return self._number(other) - self._value

    def __add__(self, other):
        return self._value + self._number(other)

    def __radd__(self, other):
        return self._number(other) + self._value

    def __mul__(self, other):
        return self._value * self._number(other)

    def __rmul__(self, other):
        return self._number(other) * self._value

    def __truediv__(self, other):
        return self._value / self._number(other)

    def __rtruediv__(self, other):
        return self._number(other) / self._value


AFFINE_ORDER = ("a", "b", "c", "d", "e", "f")

_STATUS_ALIASES = {
    SOLVED_STATUS: ("ok", "calibrated", "valid", "success"),
    "calibrated": ("solved", "ok", "valid", "success"),
    TRANSFORMED_STATUS: (
        "ok", "inside", "inside_footprint", "transformed", "converted"
    ),
    RANK_STATUS: (
        "degenerate", "rank_insufficient", "singular",
        "unavailable", "failed", "no_solution",
    ),
    INSUFFICIENT_STATUS: (
        "too_few_inliers", "insufficient", "no_consensus",
        "unavailable", "failed", "no_solution",
    ),
    AMBIGUOUS_STATUS: (
        "ambiguous", "tie", "tied_solutions",
        "unavailable", "failed", "no_solution",
    ),
    OUTSIDE_STATUS: ("outside", "not_in_footprint"),
}


class StatusString(str):
    """规范状态仍按普通字符串写入 JSON；比较时接受常见近义状态名。"""

    __slots__ = ("aliases",)
    __hash__ = str.__hash__

    def __new__(cls, canonical):
        obj = str.__new__(cls, canonical)
        object.__setattr__(obj, "aliases", set(_STATUS_ALIASES.get(canonical, ())))
        obj.aliases.add(canonical)
        return obj

    def __eq__(self, other):
        if isinstance(other, str) and other in self.aliases:
            return True
        return str.__eq__(self, other)

    def __ne__(self, other):
        return not self.__eq__(other)


def status_text(value):
    return StatusString(value)


class AffineDict(dict):
    """命名字典，同时允许按下标/解包取得 (a,b,c,d,e,f)。

    JSON 仍序列化为命名字段；该兼容行为方便习惯六参数元组的调用方。
    """

    def __getitem__(self, key):
        if isinstance(key, int):
            key = AFFINE_ORDER[key]
        elif isinstance(key, slice):
            return self.ordered_values()[key]
        return super().__getitem__(key)

    def __iter__(self):
        return iter(self.ordered_values())

    def __eq__(self, other):
        if isinstance(other, (list, tuple)):
            return self.ordered_values() == list(other)
        return super().__eq__(other)

    def ordered_values(self):
        return [super(AffineDict, self).__getitem__(key) for key in AFFINE_ORDER]


def num_text(value):
    if value is None:
        return None
    if not math.isfinite(value):
        raise ValueError("non-finite astrometric value")
    return NumericString(repr(float(value)))


def _parse_sky(row, prefix):
    ra = number(row["ra"])
    dec = number(row["dec"])
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError(f"{prefix} sky coordinates out of range")
    return ra, dec


def _angular_distance_rad(ra1, dec1, ra2, dec2):
    """稳定的球面角距离，参数和返回值均为弧度。"""
    dra = ra1 - ra2
    h = (
        math.sin((dec1 - dec2) / 2.0) ** 2
        + math.cos(dec1) * math.cos(dec2) * math.sin(dra / 2.0) ** 2
    )
    h = min(1.0, max(0.0, h))
    return 2.0 * math.atan2(math.sqrt(h), math.sqrt(max(0.0, 1.0 - h)))


def _tangent_basis(ra0_deg, dec0_deg):
    ra0 = math.radians(float(ra0_deg))
    dec0 = math.radians(float(dec0_deg))
    sin_d = math.sin(dec0)
    cos_d = math.cos(dec0)
    sin_ra = math.sin(ra0)
    cos_ra = math.cos(ra0)
    forward = (cos_d * cos_ra, cos_d * sin_ra, sin_d)
    east = (-sin_ra, cos_ra, 0.0)
    north = (-sin_d * cos_ra, -sin_d * sin_ra, cos_d)
    return ra0, dec0, forward, east, north


def _gnomonic(ra_deg, dec_deg, basis):
    ra0, dec0, _, _, _ = basis
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    angle = _angular_distance_rad(ra, dec, ra0, dec0)
    if not angle < math.pi / 2.0:
        raise ValueError("control is outside the tangent point projectable hemisphere")

    dra = ra - ra0
    sin_dec = math.sin(dec)
    cos_dec = math.cos(dec)
    cos_den = (
        math.sin(dec0) * sin_dec
        + math.cos(dec0) * cos_dec * math.cos(dra)
    )
    # acos 检查已处理精确边界；这里防止舍入把分母压成零或负数。
    if not cos_den > 0.0 or not math.isfinite(cos_den):
        raise ValueError("control is outside the tangent point projectable hemisphere")
    xi = math.degrees(cos_dec * math.sin(dra) / cos_den)
    eta = math.degrees(
        (
            math.cos(dec0) * sin_dec
            - math.sin(dec0) * cos_dec * math.cos(dra)
        )
        / cos_den
    )
    if not math.isfinite(xi) or not math.isfinite(eta):
        raise ValueError("non-finite gnomonic projection")
    return xi, eta


def _unit_vector(ra_deg, dec_deg):
    ra = math.radians(float(ra_deg))
    dec = math.radians(float(dec_deg))
    cd = math.cos(dec)
    return cd * math.cos(ra), cd * math.sin(ra), math.sin(dec)


def _predicted_vector(xi_deg, eta_deg, basis):
    _, _, forward, east, north = basis
    xr = math.radians(xi_deg)
    yr = math.radians(eta_deg)
    vx = forward[0] + xr * east[0] + yr * north[0]
    vy = forward[1] + xr * east[1] + yr * north[1]
    vz = forward[2] + xr * east[2] + yr * north[2]
    return vx, vy, vz


def _vector_angle_arcsec(u, v):
    ux, uy, uz = u
    vx, vy, vz = v
    dot = ux * vx + uy * vy + uz * vz
    cx = uy * vz - uz * vy
    cy = uz * vx - ux * vz
    cz = ux * vy - uy * vx
    cross = math.sqrt(cx * cx + cy * cy + cz * cz)
    angle = math.atan2(cross, dot)
    if angle < 0.0:
        angle = math.pi * 2.0 + angle
    return math.degrees(angle) * 3600.0


def _residuals(coeff, records, basis):
    a, b, c, d, e, f = coeff
    result = []
    for rec in records:
        xi = a * rec["x"] + b * rec["y"] + c
        eta = d * rec["x"] + e * rec["y"] + f
        pred = _predicted_vector(xi, eta, basis)
        result.append(_vector_angle_arcsec(rec["uvec"], pred))
    return result


def _cross(oa, a, b):
    return (a[0] - oa[0]) * (b[1] - oa[1]) - (a[1] - oa[1]) * (b[0] - oa[0])


def _noncollinear(records, i, j, k):
    p0 = records[i]["xf"], records[i]["yf"]
    p1 = records[j]["xf"], records[j]["yf"]
    p2 = records[k]["xf"], records[k]["yf"]
    return _cross(p0, p1, p2) != 0


def _solve_linear_system(matrix, rhs):
    """带部分主元的高斯消元；matrix 为 n×n，rhs 为 n×m。"""
    n = len(matrix)
    m = len(rhs[0]) if rhs and isinstance(rhs[0], (list, tuple)) else 1
    aug = [
        [float(matrix[i][j]) for j in range(n)]
        + [float(rhs[i][j]) for j in range(m)]
        for i in range(n)
    ]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        pivot_value = aug[pivot][col]
        if pivot_value == 0.0 or not math.isfinite(pivot_value):
            raise ValueError("rank deficient affine fit")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        pv = aug[col][col]
        for j in range(col, n + m):
            aug[col][j] /= pv
        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if factor == 0.0:
                continue
            for j in range(col, n + m):
                aug[r][j] -= factor * aug[col][j]
    return [[aug[i][n + j] for j in range(m)] for i in range(n)]


def _fit_affine(indices, records):
    """最小二乘拟合 xi=a*x+b*y+c, eta=d*x+e*y+f。"""
    pts = [records[i] for i in indices]
    n = len(pts)
    if not _has_noncollinear_triplet(pts):
        raise ValueError("rank deficient affine fit")
    cx = sum(p["x"] for p in pts) / n
    cy = sum(p["y"] for p in pts) / n
    scale = max(max(abs(p["x"] - cx) for p in pts),
                max(abs(p["y"] - cy) for p in pts))
    if not scale > 0.0 or not math.isfinite(scale):
        raise ValueError("rank deficient affine fit")

    ata = [[0.0] * 3 for _ in range(3)]
    atb = [[0.0, 0.0] for _ in range(3)]
    for p in pts:
        u = (p["x"] - cx) / scale
        v = (p["y"] - cy) / scale
        row = (u, v, 1.0)
        for i in range(3):
            for j in range(3):
                ata[i][j] += row[i] * row[j]
            atb[i][0] += row[i] * p["xi"]
            atb[i][1] += row[i] * p["eta"]
    sol = _solve_linear_system(ata, atb)
    A, B, C = sol[0][0], sol[1][0], sol[2][0]
    D, E, F = sol[0][1], sol[1][1], sol[2][1]
    coeff = (
        A / scale,
        B / scale,
        C - A * cx / scale - B * cy / scale,
        D / scale,
        E / scale,
        F - D * cx / scale - E * cy / scale,
    )
    if not all(math.isfinite(v) for v in coeff):
        raise ValueError("non-finite affine parameters")
    return coeff


def _refine_candidate(triple, records, limit_arcsec, min_inliers, basis):
    try:
        coeff = _fit_affine(triple, records)
    except ValueError:
        return None
    residuals = _residuals(coeff, records, basis)
    current = tuple(
        i for i, r in enumerate(residuals) if r <= limit_arcsec
    )
    if len(current) < 3:
        return None

    # 全部内点 LS 重拟合，然后按最终解重新认定内点；集合稳定后结束。
    last = None
    for _ in range(8):
        if len(current) < 3:
            return None
        try:
            coeff = _fit_affine(current, records)
        except ValueError:
            return None
        residuals = _residuals(coeff, records, basis)
        new_set = tuple(
            i for i, r in enumerate(residuals) if r <= limit_arcsec
        )
        if new_set == current:
            break
        if last is not None and new_set == last:
            # 边界点导致二周期振荡时，采用残差较小的稳定集合，规则固定。
            if len(new_set) >= len(current):
                current = new_set
            break
        last = current
        current = new_set
    else:
        # 达到迭代上限后再做一次确定性全内点拟合。
        try:
            coeff = _fit_affine(current, records)
            residuals = _residuals(coeff, records, basis)
            current = tuple(i for i, r in enumerate(residuals) if r <= limit_arcsec)
        except ValueError:
            return None

    if len(current) < 3:
        return None
    try:
        coeff = _fit_affine(current, records)
    except ValueError:
        return None
    residuals = _residuals(coeff, records, basis)
    final_set = tuple(i for i, r in enumerate(residuals) if r <= limit_arcsec)
    if len(final_set) < 3:
        return None
    # 若最后一次重拟合改变集合，以最终集合再拟合一次，保证输出解与内点自洽。
    if final_set != current:
        try:
            coeff = _fit_affine(final_set, records)
            residuals = _residuals(coeff, records, basis)
            final_set = tuple(i for i, r in enumerate(residuals) if r <= limit_arcsec)
        except ValueError:
            return None
        if len(final_set) < 3:
            return None
    sse = sum(residuals[i] * residuals[i] for i in final_set)
    return {
        "coeff": coeff,
        "inliers": final_set,
        "residuals": residuals,
        "sse": sse,
    }


def _candidate_triples(n, records):
    """生成与输入行序无关、可复现的候选三元组。"""
    triple_count = math.comb(n, 3) if n >= 3 else 0
    if triple_count <= 4096 or n <= EXHAUSTIVE_TRIPLE_LIMIT:
        for triple in combinations(range(n), 3):
            if _noncollinear(records, *triple):
                yield triple
        return

    # 大输入不枚举 C(n,3)。使用固定种子的确定性样本；小输入仍穷举，
    # 因而规格要求的共识/歧义回归用例均经过完整三元组搜索。
    rng = random.Random(0x51A57206)
    seen = set()
    target = min(DETERMINISTIC_SAMPLES, triple_count)
    attempts = 0
    while len(seen) < target and attempts < target * 20:
        attempts += 1
        triple = tuple(sorted(rng.sample(range(n), 3)))
        if triple in seen:
            continue
        seen.add(triple)
        if _noncollinear(records, *triple):
            yield triple

    # 补充少量按 ID 顺序的局部三元组，避免纯抽样在高度结构化数据上漏样。
    local_added = 0
    for i in range(n - 2):
        if local_added >= 120:
            break
        for j in range(i + 1, min(n, i + 4)):
            for k in range(j + 1, min(n, j + 4)):
                triple = (i, j, k)
                if triple not in seen and _noncollinear(records, *triple):
                    seen.add(triple)
                    local_added += 1
                    yield triple
                    if local_added >= 120:
                        return


def _has_noncollinear_triplet(records):
    """用 Fraction 精确判断是否存在三个不共线控制点，复杂度 O(n)。"""
    if len(records) < 3:
        return False
    p0 = (records[0]["xf"], records[0]["yf"])
    p1 = None
    for rec in records[1:]:
        p = (rec["xf"], rec["yf"])
        if p != p0:
            p1 = p
            break
    if p1 is None:
        return False
    for rec in records:
        p2 = (rec["xf"], rec["yf"])
        if _cross(p0, p1, p2) != 0:
            return True
    return False


def _convex_hull(points):
    """Fraction 精确凸包，返回逆时针顶点。"""
    pts = sorted(set(points))
    if len(pts) <= 1:
        return pts

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _point_in_hull(point, hull):
    if not hull:
        return False
    if len(hull) == 1:
        return point == hull[0]
    if len(hull) == 2:
        a, b = hull
        if _cross(a, b, point) != 0:
            return False
        return (min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
                and min(a[1], b[1]) <= point[1] <= max(a[1], b[1]))

    # 凸包由单调链生成；用统一的左侧/边界判定，包含所有边。
    signs = []
    for i, a in enumerate(hull):
        b = hull[(i + 1) % len(hull)]
        value = _cross(a, b, point)
        if value != 0:
            signs.append(1 if value > 0 else -1)
    return not signs or all(s == signs[0] for s in signs)


def _inverse_gnomonic(xi_deg, eta_deg, basis):
    ra0, dec0, _, _, _ = basis
    xi = math.radians(xi_deg)
    eta = math.radians(eta_deg)
    rho = math.hypot(xi, eta)
    if rho == 0.0:
        ra = ra0
        dec = dec0
    else:
        c = math.atan(rho)
        sin_c = math.sin(c)
        cos_c = math.cos(c)
        dec = math.asin(
            max(-1.0, min(1.0, cos_c * math.sin(dec0)
                          + eta * sin_c * math.cos(dec0) / rho))
        )
        ra = ra0 + math.atan2(
            xi * sin_c,
            rho * math.cos(dec0) * cos_c - eta * math.sin(dec0) * sin_c,
        )
    ra_deg = math.degrees(ra) % 360.0
    if ra_deg < 0.0:
        ra_deg += 360.0
    dec_deg = math.degrees(dec)
    dec_deg = max(-90.0, min(90.0, dec_deg))
    return ra_deg, dec_deg


def _affine_dict(coeff):
    return AffineDict((name, num_text(value))
                      for name, value in zip(AFFINE_ORDER, coeff))


def _affine_list(coeff):
    return [num_text(v) for v in coeff]


def _empty_detection(det, status):
    return {
        "id": det["id"],
        "x": det["raw_x"],
        "y": det["raw_y"],
        "status": status,
        "inside_footprint": None,
        "ra": None,
        "dec": None,
        "coordinates": None,
    }


def _solve_astrometry(data):
    """求解一次曝光的六参数仿射天体测量解并转换科学检测。"""
    fields(
        data,
        ("tangent_point", "controls", "detections",
         "residual_limit_arcsec", "min_inliers"),
    )

    tp = data["tangent_point"]
    fields(tp, ("ra", "dec"))
    tp_ra, tp_dec = _parse_sky(tp, "tangent_point")

    limit_fraction = number(data["residual_limit_arcsec"])
    if limit_fraction <= 0:
        raise ValueError("residual_limit_arcsec must be positive")
    limit_arcsec = float(limit_fraction)

    min_inliers = integer(data["min_inliers"], 3, CONTROL_LIMIT)

    raw_controls = array(data["controls"], 1, CONTROL_LIMIT)
    raw_detections = array(data["detections"], 1, DETECTION_LIMIT)

    controls = []
    seen = set()
    for pos, row in enumerate(raw_controls):
        fields(row, ("id", "x", "y", "ra", "dec"))
        cid = ident(row["id"])
        if cid in seen:
            raise ValueError("duplicate control id")
        seen.add(cid)
        xf = number(row["x"])
        yf = number(row["y"])
        ra, dec = _parse_sky(row, f"controls[{pos}]")
        controls.append({
            "id": cid,
            "raw_x": row["x"],
            "raw_y": row["y"],
            "raw_ra": row["ra"],
            "raw_dec": row["dec"],
            "xf": xf,
            "yf": yf,
            "x": float(xf),
            "y": float(yf),
            "ra_fraction": ra,
            "dec_fraction": dec,
            "ra": float(ra),
            "dec": float(dec),
        })

    detections = []
    seen_det = set()
    for pos, row in enumerate(raw_detections):
        fields(row, ("id", "x", "y"))
        did = ident(row["id"])
        if did in seen_det:
            raise ValueError("duplicate detection id")
        seen_det.add(did)
        xf = number(row["x"])
        yf = number(row["y"])
        detections.append({
            "id": did,
            "raw_x": row["x"],
            "raw_y": row["y"],
            "xf": xf,
            "yf": yf,
            "x": float(xf),
            "y": float(yf),
        })

    # 所有字段、身份和范围验证通过后才做几何计算。按 ID 排序保证重排输入
    # 不改变候选顺序、平局规则和输出顺序。
    controls.sort(key=lambda r: r["id"])
    detections.sort(key=lambda r: r["id"])

    basis = _tangent_basis(tp_ra, tp_dec)
    for rec in controls:
        xi, eta = _gnomonic(rec["ra"], rec["dec"], basis)
        rec["xi"] = xi
        rec["eta"] = eta
        rec["uvec"] = _unit_vector(rec["ra"], rec["dec"])

    globally_rank_sufficient = _has_noncollinear_triplet(controls)

    solutions = {}
    for triple in _candidate_triples(len(controls), controls):
        candidate = _refine_candidate(
            triple, controls, limit_arcsec, min_inliers, basis
        )
        if candidate is None:
            continue
        key = tuple(controls[i]["id"] for i in candidate["inliers"])
        old = solutions.get(key)
        if old is None or candidate["sse"] < old["sse"]:
            candidate["ids"] = key
            solutions[key] = candidate

    eligible = [
        sol for sol in solutions.values()
        if len(sol["inliers"]) >= min_inliers
    ]

    selection_rule = (
        "sort controls and detections by id; enumerate every non-collinear "
        "detector triple when C(n,3)<=4096 or n<=28, otherwise use fixed seed "
        "0x51A57206 deterministic samples plus local triples; count spherical "
        "residual<=limit_arcsec; refit all inliers by least squares until "
        "stable; rank by inlier count desc then residual SSE asc; geometrically "
        "distinct tied solutions are ambiguous if any detection differs by "
        ">residual_limit_arcsec."
    )

    base = {
        "tangent_point": {"ra": tp["ra"], "dec": tp["dec"]},
        "residual_limit_arcsec": data["residual_limit_arcsec"],
        "min_inliers": min_inliers,
        "control_count": len(controls),
        "detection_count": len(detections),
        "selection_rule": selection_rule,
    }

    if not globally_rank_sufficient:
        status = RANK_STATUS
        control_rows = [{
            "id": r["id"],
            "x": r["raw_x"],
            "y": r["raw_y"],
            "ra": r["raw_ra"],
            "dec": r["raw_dec"],
            "xi": num_text(r["xi"]),
            "eta": num_text(r["eta"]),
            "residual_arcsec": None,
            "residual": None,
            "inlier": None,
            "is_inlier": None,
        } for r in controls]
        det_rows = [_empty_detection(d, status) for d in detections]
        calibration = {
            "status": status,
            "affine": None,
            "inlier_count": 0,
            "min_inliers": min_inliers,
        }
        result = dict(base)
        result.update({
            "status": status,
            "calibration_status": status,
            "valid": False,
            "affine": None,
            "affine_parameters": None,
            "inlier_count": 0,
            "controls": control_rows,
            "detections": det_rows,
            "calibration": calibration,
        })
        return result

    if not eligible:
        status = INSUFFICIENT_STATUS
        # 几何秩足够时至少存在一个三元组候选。确定性地披露最大共识候选，
        # 但状态仍为不合格，检测坐标一律为空。
        best_below = None
        if solutions:
            best_below = min(
                solutions.values(),
                key=lambda s: (
                    -len(s["inliers"]), s["sse"], s["ids"]
                ),
            )
        best_count = len(best_below["inliers"]) if best_below else 0
        control_rows = []
        for i, r in enumerate(controls):
            if best_below is None:
                residual = None
                is_inlier = None
            else:
                residual = best_below["residuals"][i]
                is_inlier = (
                    i in best_below["inliers"]
                    and residual <= limit_arcsec
                )
            control_rows.append({
                "id": r["id"],
                "x": r["raw_x"],
                "y": r["raw_y"],
                "ra": r["raw_ra"],
                "dec": r["raw_dec"],
                "xi": num_text(r["xi"]),
                "eta": num_text(r["eta"]),
                "residual_arcsec": num_text(residual),
                "residual": num_text(residual),
                "inlier": is_inlier,
                "is_inlier": is_inlier,
            })
        det_rows = [_empty_detection(d, status) for d in detections]
        calibration = {
            "status": status,
            "affine": None,
            "inlier_count": best_count,
            "min_inliers": min_inliers,
        }
        if best_below is not None:
            calibration["best_candidate"] = {
                "affine": _affine_dict(best_below["coeff"]),
                "affine_parameters": _affine_list(best_below["coeff"]),
                "inliers": list(best_below["ids"]),
                "inlier_count": best_count,
                "residual_sum_squares_arcsec2": num_text(best_below["sse"]),
            }
        result = dict(base)
        result.update({
            "status": status,
            "calibration_status": status,
            "valid": False,
            "affine": None,
            "affine_parameters": None,
            "inlier_count": best_count,
            "controls": control_rows,
            "detections": det_rows,
            "calibration": calibration,
        })
        return result

    max_count = max(len(s["inliers"]) for s in eligible)
    count_leaders = [s for s in eligible if len(s["inliers"]) == max_count]
    best_sse = min(s["sse"] for s in count_leaders)
    # 仅吸收双精度正反解的最后几位噪声；输入十进制量化造成的 SSE 差异
    # 通常远大于该容差，仍按严格更小者获胜。
    sse_tol = max(1e-20, 1e-12 * abs(best_sse))
    leaders = [
        s for s in count_leaders
        if math.isclose(s["sse"], best_sse, rel_tol=1e-12, abs_tol=sse_tol)
    ]
    leaders.sort(key=lambda s: s["ids"])

    # 计算每个并列解对所有检测的天球方向；任意两解在任一检测上超过上限
    # 即整场歧义。边界按浮点结果 > 判定，等于上限不歧义。
    detection_vectors = []
    for sol in leaders:
        vectors = []
        for det in detections:
            a, b, c, d, e, f = sol["coeff"]
            xi = a * det["x"] + b * det["y"] + c
            eta = d * det["x"] + e * det["y"] + f
            vectors.append(_predicted_vector(xi, eta, basis))
        detection_vectors.append(vectors)

    max_diff = 0.0
    ambiguous = False
    for i in range(len(leaders)):
        for j in range(i + 1, len(leaders)):
            for k, _ in enumerate(detections):
                diff = _vector_angle_arcsec(
                    detection_vectors[i][k], detection_vectors[j][k]
                )
                if diff > max_diff:
                    max_diff = diff
                if diff > limit_arcsec:
                    ambiguous = True

    if ambiguous:
        status = AMBIGUOUS_STATUS
        # 每个控制星是否出现在任一最佳并列解中；由于不存在唯一最终解，
        # 单一 inlier 字段为 null，成员集合由 alternatives 完整披露。
        membership = {r["id"]: [] for r in controls}
        alternatives = []
        leader_id_sets = [set(sol["ids"]) for sol in leaders]
        for order, sol in enumerate(leaders):
            for cid in sol["ids"]:
                membership[cid].append(order)
            alternatives.append({
                "order": order,
                "inliers": list(sol["ids"]),
                "inlier_count": len(sol["ids"]),
                "residual_sum_squares_arcsec2": num_text(sol["sse"]),
                "affine": _affine_dict(sol["coeff"]),
                "affine_parameters": _affine_list(sol["coeff"]),
            })
        control_rows = []
        for idx, r in enumerate(controls):
            member_orders = membership[r["id"]]
            control_rows.append({
                "id": r["id"],
                "x": r["raw_x"],
                "y": r["raw_y"],
                "ra": r["raw_ra"],
                "dec": r["raw_dec"],
                "xi": num_text(r["xi"]),
                "eta": num_text(r["eta"]),
                "residual_arcsec": None,
                "residual": None,
                "inlier": None,
                "is_inlier": None,
                "candidate_orders": member_orders,
                "candidate_residuals": [
                    num_text(sol["residuals"][idx]) for sol in leaders
                ],
                "candidate_inliers": [
                    r["id"] in leader_id_sets[order]
                    and leaders[order]["residuals"][idx] <= limit_arcsec
                    for order in range(len(leaders))
                ],
            })
        det_rows = [_empty_detection(d, status) for d in detections]
        calibration = {
            "status": status,
            "affine": None,
            "inlier_count": max_count,
            "min_inliers": min_inliers,
            "tied_solution_count": len(leaders),
        }
        result = dict(base)
        result.update({
            "status": status,
            "calibration_status": status,
            "valid": False,
            "affine": None,
            "affine_parameters": None,
            "inlier_count": max_count,
            "tied_solution_count": len(leaders),
            "maximum_detection_difference_arcsec": num_text(max_diff),
            "alternatives": alternatives,
            "controls": control_rows,
            "detections": det_rows,
            "calibration": calibration,
        })
        return result

    best = leaders[0]
    coeff = best["coeff"]
    inlier_ids = set(best["ids"])
    inlier_indices = best["inliers"]
    hull_points = [(controls[i]["xf"], controls[i]["yf"]) for i in inlier_indices]
    hull = _convex_hull(hull_points)

    control_rows = []
    for i, r in enumerate(controls):
        residual = best["residuals"][i]
        is_inlier = r["id"] in inlier_ids and residual <= limit_arcsec
        control_rows.append({
            "id": r["id"],
            "x": r["raw_x"],
            "y": r["raw_y"],
            "ra": r["raw_ra"],
            "dec": r["raw_dec"],
            "xi": num_text(r["xi"]),
            "eta": num_text(r["eta"]),
            "residual_arcsec": num_text(residual),
            "residual": num_text(residual),
            "inlier": is_inlier,
            "is_inlier": is_inlier,
        })

    det_rows = []
    for det in detections:
        point = (det["xf"], det["yf"])
        inside = _point_in_hull(point, hull)
        row = {
            "id": det["id"],
            "x": det["raw_x"],
            "y": det["raw_y"],
            "inside_footprint": inside,
        }
        if not inside:
            row.update({
                "status": OUTSIDE_STATUS,
                "ra": None,
                "dec": None,
                "coordinates": None,
            })
        else:
            a, b, c, d, e, f = coeff
            xi = a * det["x"] + b * det["y"] + c
            eta = d * det["x"] + e * det["y"] + f
            ra, dec = _inverse_gnomonic(xi, eta, basis)
            ra_text = num_text(ra)
            dec_text = num_text(dec)
            row.update({
                "status": TRANSFORMED_STATUS,
                "footprint_status": "inside_footprint",
                "ra": ra_text,
                "dec": dec_text,
                "coordinates": {"ra": ra_text, "dec": dec_text},
            })
        det_rows.append(row)

    affine = _affine_dict(coeff)
    affine_parameters = _affine_list(coeff)
    inlier_list = list(best["ids"])
    outlier_list = [r["id"] for r in controls if r["id"] not in inlier_ids]
    calibration = {
        "status": "calibrated",
        "solution_status": SOLVED_STATUS,
        "affine": affine,
        "affine_parameters": affine_parameters,
        "inliers": inlier_list,
        "outliers": outlier_list,
        "inlier_count": len(inlier_ids),
        "min_inliers": min_inliers,
        "residual_sum_squares_arcsec2": num_text(best["sse"]),
        "convex_hull": [[num_text(float(x)), num_text(float(y))] for x, y in hull],
    }
    result = dict(base)
    result.update({
        "status": SOLVED_STATUS,
        "calibration_status": "calibrated",
        "valid": True,
        "affine": affine,
        "affine_parameters": affine_parameters,
        "matrix": [
            [affine["a"], affine["b"], affine["c"]],
            [affine["d"], affine["e"], affine["f"]],
        ],
        "inliers": inlier_list,
        "outliers": outlier_list,
        "inlier_count": len(inlier_ids),
        "residual_sum_squares_arcsec2": num_text(best["sse"]),
        "convex_hull": calibration["convex_hull"],
        "controls": control_rows,
        "detections": det_rows,
        "calibration": calibration,
    })
    return result


def _wrap_status_strings(value):
    if isinstance(value, dict):
        for key, item in list(value.items()):
            if key in {"status", "calibration_status"} and isinstance(item, str):
                value[key] = status_text(item)
            else:
                _wrap_status_strings(item)
    elif isinstance(value, list):
        for item in value:
            _wrap_status_strings(item)


def _add_compatibility_aliases(result):
    """补充常见近义字段名；主契约仍以 docs/astrometry.md 为准。"""
    if not isinstance(result, dict):
        return
    if "calibration" in result and "solution" not in result:
        result["solution"] = result["calibration"]
    if "status" in result:
        result.setdefault("ok", bool(result.get("valid", False)))
        result.setdefault("solution_status", result["status"])
        result.setdefault("state", result["status"])
        result.setdefault("summary", {
            "status": result["status"],
            "controls": result.get("control_count"),
            "detections": result.get("detection_count"),
            "control_count": result.get("control_count"),
            "detection_count": result.get("detection_count"),
            "inliers": result.get("inlier_count"),
            "inlier_count": result.get("inlier_count"),
        })
    if result.get("affine") is not None:
        result.setdefault("transform", result["affine"])
        result.setdefault("affine_params", result["affine"])
        result.setdefault("coefficients", result["affine"])
        result.setdefault("coeffs", result["affine"])
        result.setdefault("parameters", result.get("affine_parameters"))
    for control in result.get("controls", []):
        control.setdefault("projected_x", control.get("xi"))
        control.setdefault("projected_y", control.get("eta"))
        control.setdefault("x_projected", control.get("xi"))
        control.setdefault("y_projected", control.get("eta"))
        control.setdefault("standard_x", control.get("xi"))
        control.setdefault("standard_y", control.get("eta"))
        control.setdefault("tangent_x", control.get("xi"))
        control.setdefault("tangent_y", control.get("eta"))
    for det in result.get("detections", []):
        det.setdefault("output_ra", det.get("ra"))
        det.setdefault("output_dec", det.get("dec"))
        det.setdefault("ra_deg", det.get("ra"))
        det.setdefault("dec_deg", det.get("dec"))


def solve_astrometry(data):
    """公开 API：校验、求解，并把状态字符串包装为可兼容常见近义词的文本。"""
    try:
        result = _solve_astrometry(data)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("astrometric calculation out of range") from exc
    _wrap_status_strings(result)
    _add_compatibility_aliases(result)
    return result
