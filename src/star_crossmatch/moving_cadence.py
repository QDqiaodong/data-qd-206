"""三帧巡天的移动轨迹审计与静态测光交接.

轨迹候选必须由早、中、晚三帧各一条检测组成。速度闭区间和中间点线性插值
容差全部用 ``Fraction`` 精确判定;代价是中间预测点到检测点的欧氏距离平方,
序列化为约分有理数字符串。候选按共享检测形成冲突分量,分量内参与候选的
不同检测超过 15 条时整批拒绝。

静态测光不是给原结果补标签:每个全批最佳轨迹方案都使用的检测先从静态输
入剔除;仅部分最佳方案使用的检测保留,但其静态配对不得作为零点锚点。其余
流程完全复用 ``process_exposures`` 的球面候选、整数角距代价、一对一几何最
优、锚点/零点和星等破平局规则。
"""
from fractions import Fraction
from functools import lru_cache

from .checks import fields, array, ident, number
from .exposures import (
    _decimal_str,
    _process_exposure,
    _validate,
    process_exposures,
)

DECIMAL_SCALE = 1_000_000
MAX_COMPONENT_DETECTIONS = 15

_STATE_LABEL = {
    "always": "always_used",
    "partial": "partially_used",
    "never": "never_used",
}


class _DSU:
    """按检测身份合并候选冲突分量的轻量并查集。"""

    def __init__(self, n):
        self.parent = list(range(n))
        self.size = [1] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return ra
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]
        return ra


def _micro(value):
    """现有十进制语法最多 6 位小数,可无损转成百万分位整数。"""
    return int(value * DECIMAL_SCALE)


def rational_string(value):
    """把精确分数写成约分有理数字符串;整数不输出分母 1。"""
    if not isinstance(value, Fraction):
        value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _identity_ref(key):
    exposure, observation = key
    return {
        "exposure": exposure,
        "observation": observation,
        "identity": [exposure, observation],
    }


def _better_plan(candidate, current):
    """轨迹数更多优先,总代价更小优先,最后取稳定的最小见证键。"""
    if current is None:
        return True
    if candidate[0] != current[0]:
        return candidate[0] > current[0]
    if candidate[1] != current[1]:
        return candidate[1] < current[1]
    return candidate[2] < current[2]


def _same_objective(left, right):
    return left[0] == right[0] and left[1] == right[1]


class _ComponentSolver:
    """不超过 15 个检测位的分量集合打包 DP。

    ``best(available)`` 返回 ``(轨迹数, 总代价, 规范化候选下标见证)``。
    状态是仍可使用的检测位掩码;选择一条候选后删除它占用的全部检测位,因
    此同一候选不可能被再次使用,见证内也不可能复用检测。
    """

    def __init__(self, candidates, size):
        self.candidates = candidates
        self.size = size
        self.full = (1 << size) - 1

    @lru_cache(maxsize=None)
    def best(self, available):
        best = (0, Fraction(0), ())
        for index, candidate in enumerate(self.candidates):
            mask = candidate["mask"]
            if mask & ~available:
                continue
            count, cost, witness = self.best(available ^ mask)
            plan = (
                count + 1,
                cost + candidate["cost"],
                tuple(sorted((index,) + witness)),
            )
            if _better_plan(plan, best):
                best = plan
        return best

    def required_detection(self, bit):
        best = None
        for index, candidate in enumerate(self.candidates):
            if not candidate["mask"] & bit:
                continue
            count, cost, witness = self.best(self.full ^ candidate["mask"])
            plan = (
                count + 1,
                cost + candidate["cost"],
                tuple(sorted((index,) + witness)),
            )
            if _better_plan(plan, best):
                best = plan
        return best

    def forbidden_detection(self, bit):
        return self.best(self.full ^ bit)


def _validate_input(data):
    fields(
        data,
        (
            "photometry",
            "epochs",
            "positions",
            "min_speed",
            "max_speed",
            "middle_tolerance",
        ),
    )

    min_speed = number(data["min_speed"])
    max_speed = number(data["max_speed"])
    tolerance = number(data["middle_tolerance"])
    if min_speed <= 0 or max_speed <= 0 or min_speed > max_speed:
        raise ValueError("requires 0 < min_speed <= max_speed")
    if tolerance < 0:
        raise ValueError("middle_tolerance must be non-negative")

    threshold, anchor_tolerance, reference, exposures = _validate(
        data["photometry"]
    )
    if len(exposures) != 3:
        raise ValueError("reduce_moving_cadence requires exactly 3 exposures")
    if not 1 <= len(reference) <= 256:
        raise ValueError("reference count 1..256 required")
    for exposure in exposures:
        if not 1 <= len(exposure["observations"]) <= 120:
            raise ValueError("each exposure requires 1..120 detections")

    epoch_rows = array(data["epochs"], 3, 3)
    epochs = {}
    for row in epoch_rows:
        fields(row, ("exposure", "mjd"))
        exposure_id = ident(row["exposure"])
        if exposure_id in epochs:
            raise ValueError("duplicate epoch exposure")
        mjd = number(row["mjd"])
        if mjd < 0:
            raise ValueError("mjd must be non-negative")
        epochs[exposure_id] = mjd

    exposure_ids = {exposure["id"] for exposure in exposures}
    if set(epochs) != exposure_ids:
        raise ValueError("epochs must contain exactly one row per exposure")
    if len(set(epochs.values())) != 3:
        raise ValueError("mjd values must be distinct")

    observation_keys = {
        (exposure["id"], observation["id"])
        for exposure in exposures
        for observation in exposure["observations"]
    }
    total = len(observation_keys)
    position_rows = array(data["positions"], total, total)
    positions = {}
    for row in position_rows:
        fields(row, ("exposure", "observation", "x", "y"))
        key = (ident(row["exposure"]), ident(row["observation"]))
        if key not in observation_keys:
            raise ValueError("position refers to missing detection")
        if key in positions:
            raise ValueError("duplicate position")
        x = number(row["x"])
        y = number(row["y"])
        positions[key] = (x, y, _micro(x), _micro(y))

    if set(positions) != observation_keys:
        raise ValueError("every detection requires exactly one position")

    return (
        threshold,
        anchor_tolerance,
        reference,
        exposures,
        min_speed,
        max_speed,
        tolerance,
        epochs,
        positions,
    )


def _build_candidates(exposures, epochs, positions, min_speed, max_speed, tolerance):
    """生成精确入围候选,并在候选冲突分量超过 15 个检测时立即拒绝。"""
    ordered = sorted(exposures, key=lambda exposure: epochs[exposure["id"]])
    early_id, middle_id, late_id = (exposure["id"] for exposure in ordered)
    epoch_micro = {
        exposure_id: _micro(mjd) for exposure_id, mjd in epochs.items()
    }

    all_keys = sorted(positions)
    gid = {key: index for index, key in enumerate(all_keys)}
    role_keys = {
        "early": sorted(key for key in all_keys if key[0] == early_id),
        "middle": sorted(key for key in all_keys if key[0] == middle_id),
        "late": sorted(key for key in all_keys if key[0] == late_id),
    }

    dsu = _DSU(len(all_keys))
    min_speed_micro = _micro(min_speed)
    max_speed_micro = _micro(max_speed)
    tolerance_micro = _micro(tolerance)
    scale_squared = DECIMAL_SCALE * DECIMAL_SCALE
    raw_candidates = []

    for early_key in role_keys["early"]:
        _, _, xe, ye = positions[early_key]
        early_time = epoch_micro[early_id]
        early_gid = gid[early_key]

        for late_key in role_keys["late"]:
            _, _, xl, yl = positions[late_key]
            late_time = epoch_micro[late_id]
            late_gid = gid[late_key]

            dx = xl - xe
            dy = yl - ye
            dt = late_time - early_time
            endpoint_squared = dx * dx + dy * dy

            # speed = sqrt(endpoint_squared) / dt;速度边界乘 1e6 后得到
            # 这里的 1e12 * endpoint_squared 比较,全程无浮点。
            if scale_squared * endpoint_squared < (
                min_speed_micro * min_speed_micro * dt * dt
            ):
                continue
            if scale_squared * endpoint_squared > (
                max_speed_micro * max_speed_micro * dt * dt
            ):
                continue

            for middle_key in role_keys["middle"]:
                _, _, xm, ym = positions[middle_key]
                middle_time = epoch_micro[middle_id]
                middle_gid = gid[middle_key]

                dt_middle = middle_time - early_time
                numerator_x = (xm - xe) * dt - dx * dt_middle
                numerator_y = (ym - ye) * dt - dy * dt_middle
                residual_squared = (
                    numerator_x * numerator_x + numerator_y * numerator_y
                )
                if residual_squared > (
                    tolerance_micro * tolerance_micro * dt * dt
                ):
                    continue

                cost = Fraction(
                    residual_squared,
                    dt * dt * scale_squared,
                )
                raw_candidates.append(
                    {
                        "keys": (early_key, middle_key, late_key),
                        "gids": (early_gid, middle_gid, late_gid),
                        "cost": cost,
                        "speed_squared": Fraction(endpoint_squared, dt * dt),
                    }
                )
                dsu.union(early_gid, middle_gid)
                dsu.union(middle_gid, late_gid)
                root = dsu.find(early_gid)
                if dsu.size[root] > MAX_COMPONENT_DETECTIONS:
                    raise ValueError(
                        "moving conflict component contains more than 15 detections"
                    )

    raw_candidates.sort(key=lambda candidate: candidate["keys"])
    return ordered, all_keys, gid, dsu, raw_candidates


def reduce_moving_cadence(data):
    """归约三帧连续巡天,返回轨迹审计和可交接恒星组的静态测光结果。"""
    (
        threshold,
        anchor_tolerance,
        reference,
        exposures,
        min_speed,
        max_speed,
        tolerance,
        epochs,
        positions,
    ) = _validate_input(data)

    ordered, all_keys, gid, dsu, raw_candidates = _build_candidates(
        exposures, epochs, positions, min_speed, max_speed, tolerance
    )

    grouped = {}
    for candidate in raw_candidates:
        root = dsu.find(candidate["gids"][0])
        grouped.setdefault(root, []).append(candidate)

    component_builders = []
    for candidates in grouped.values():
        member_gids = sorted(
            {member for candidate in candidates for member in candidate["gids"]}
        )
        if len(member_gids) > MAX_COMPONENT_DETECTIONS:
            raise ValueError(
                "moving conflict component contains more than 15 detections"
            )
        local_index = {global_id: index for index, global_id in enumerate(member_gids)}
        local_candidates = []
        for index, candidate in enumerate(candidates):
            mask = 0
            for global_id in candidate["gids"]:
                mask |= 1 << local_index[global_id]
            local_candidates.append({**candidate, "mask": mask, "local_index": index})
        component_builders.append(
            {
                "member_gids": member_gids,
                "candidates": local_candidates,
            }
        )

    component_builders.sort(
        key=lambda component: tuple(
            all_keys[global_id] for global_id in component["member_gids"]
        )
    )

    components = []
    candidate_records = []
    records_by_gid = {}
    always_keys = set()
    partial_keys = set()
    batch_count = 0
    batch_cost = Fraction(0)

    for component_number, builder in enumerate(component_builders, start=1):
        component_id = f"C{component_number}"
        member_gids = builder["member_gids"]
        local_candidates = builder["candidates"]
        size = len(member_gids)
        solver = _ComponentSolver(local_candidates, size)
        unconstrained = solver.best(solver.full)
        batch_count += unconstrained[0]
        batch_cost += unconstrained[1]

        def plan_payload(plan):
            if plan is None:
                return None
            count, total_cost, witness = plan
            witness_triples = []
            identity_witness = []
            for candidate_index in witness:
                candidate = local_candidates[candidate_index]
                identities = [
                    _identity_ref(all_keys[global_id])
                    for global_id in candidate["gids"]
                ]
                witness_triples.append(identities)
                identity_witness.append(
                    [identity["identity"] for identity in identities]
                )
            cost_text = rational_string(total_cost)
            return {
                "feasible": True,
                "trajectory_count": count,
                "count": count,
                "total_cost": cost_text,
                "cost": cost_text,
                "witness": witness_triples,
                "identity_witness": identity_witness,
            }

        best_payload = plan_payload(unconstrained)
        detection_refs = [_identity_ref(all_keys[g]) for g in member_gids]
        component_candidate_records = []

        for local_index_value, candidate in enumerate(local_candidates):
            candidate_id = f"{component_id}-T{local_index_value + 1}"
            identities = [
                _identity_ref(all_keys[global_id])
                for global_id in candidate["gids"]
            ]
            record = {
                "id": candidate_id,
                "candidate_id": candidate_id,
                "component": component_id,
                "early": identities[0],
                "middle": identities[1],
                "late": identities[2],
                "cost": rational_string(candidate["cost"]),
                "speed_squared": rational_string(candidate["speed_squared"]),
            }
            component_candidate_records.append(record)
            candidate_records.append(record)

        component_record = {
            "id": component_id,
            "component": component_id,
            "component_id": component_id,
            "component_index": component_number,
            "detection_count": size,
            "detections": detection_refs,
            "candidate_count": len(local_candidates),
            "candidates": component_candidate_records,
            "unconstrained_best": best_payload,
            "best": best_payload,
        }
        components.append(component_record)

        for local_detection, global_id in enumerate(member_gids):
            bit = 1 << local_detection
            required_plan = solver.required_detection(bit)
            forbidden_plan = solver.forbidden_detection(bit)
            required_payload = plan_payload(required_plan)
            forbidden_payload = plan_payload(forbidden_plan)

            can_require_best = _same_objective(required_plan, unconstrained)
            can_forbid_best = _same_objective(forbidden_plan, unconstrained)
            if can_require_best and can_forbid_best:
                state = "partial"
                partial_keys.add(all_keys[global_id])
            elif can_require_best:
                state = "always"
                always_keys.add(all_keys[global_id])
            elif can_forbid_best:
                state = "never"
            else:  # pragma: no cover - 数学上不可能,保留以尽早暴露求解错误
                raise RuntimeError("inconsistent moving-cadence optimal states")

            key = all_keys[global_id]
            identity = _identity_ref(key)
            if state == "always":
                static_action = "remove"
            elif state == "partial":
                static_action = "retain_but_block_anchor"
            else:
                static_action = "retain"
            counterfactuals = {
                "require_used": required_payload,
                "forbid_used": forbidden_payload,
                "if_used": required_payload,
                "if_unused": forbidden_payload,
            }
            records_by_gid[global_id] = {
                **identity,
                "component": component_id,
                "component_id": component_id,
                "component_index": component_number,
                "in_candidates": True,
                "state": state,
                "status": _STATE_LABEL[state],
                "static_action": static_action,
                "removed_from_static": state == "always",
                "anchor_blocked": state == "partial",
                "used_by_every_best_plan": state == "always",
                "used_by_some_best_plan": state in ("always", "partial"),
                "counterfactuals": counterfactuals,
            }

    for global_id, key in enumerate(all_keys):
        if global_id in records_by_gid:
            continue
        identity = _identity_ref(key)
        records_by_gid[global_id] = {
            **identity,
            "component": None,
            "component_id": None,
            "component_index": None,
            "in_candidates": False,
            "state": "not_candidate",
            "status": "not_in_candidate",
            "static_action": "retain",
            "removed_from_static": False,
            "anchor_blocked": False,
            "used_by_every_best_plan": False,
            "used_by_some_best_plan": False,
            "counterfactuals": None,
        }

    audit = [records_by_gid[global_id] for global_id in range(len(all_keys))]

    # 没有 always/partial 限制时,静态结果必须与原 API 完全一致;直接调用原
    # 入口也避免复制任何边界行为。有限制时从规范化输入构造剔除后的帧。
    if not always_keys and not partial_keys:
        static_result = process_exposures(data["photometry"])
    else:
        static_exposures = []
        for exposure in sorted(exposures, key=lambda item: item["id"]):
            blocked = {
                observation_id
                for exposure_id, observation_id in partial_keys
                if exposure_id == exposure["id"]
            }
            observations = [
                observation
                for observation in exposure["observations"]
                if (exposure["id"], observation["id"]) not in always_keys
            ]
            static_exposures.append(
                {"id": exposure["id"], "observations": observations}
            )
        static_outputs = [
            _process_exposure(
                exposure,
                threshold,
                anchor_tolerance,
                reference,
                anchor_blocked=blocked,
            )
            for exposure in static_exposures
        ]
        static_result = {"exposures": static_outputs}

    removed_refs = [_identity_ref(key) for key in sorted(always_keys)]
    blocked_refs = [_identity_ref(key) for key in sorted(partial_keys)]
    retained_refs = [
        _identity_ref(key) for key in all_keys if key not in always_keys
    ]

    frames = [
        {
            "role": role,
            "exposure": exposure["id"],
            "mjd": _decimal_str(epochs[exposure["id"]]),
        }
        for role, exposure in zip(("early", "middle", "late"), ordered)
    ]
    global_cost = rational_string(batch_cost)
    global_best = {
        "trajectory_count": batch_count,
        "count": batch_count,
        "total_cost": global_cost,
        "cost": global_cost,
    }
    parameters = {
        "min_speed": _decimal_str(min_speed),
        "max_speed": _decimal_str(max_speed),
        "middle_tolerance": _decimal_str(tolerance),
    }
    handoff = {
        "removed_detections": removed_refs,
        "anchor_blocked_detections": blocked_refs,
        "retained_detections": retained_refs,
    }
    return {
        "schema": "star_crossmatch.moving_cadence/v1",
        "api": "reduce_moving_cadence",
        "frames": frames,
        "epochs": frames,
        "parameters": parameters,
        "global_best": global_best,
        "components": components,
        "candidates": candidate_records,
        "detections": audit,
        "audit": audit,
        "exposures": static_result["exposures"],
        "static_photometry": static_result,
        "photometry": static_result,
        "static_handoff": handoff,
    }
