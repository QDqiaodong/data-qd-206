"""离线测光发布件与可复核变暗判据。

本模块只依赖标准库,所有星等、零点与差值均以 :class:`fractions.Fraction`
按输入十进制的精确值计算,不使用浮点近似决定相等、定标矛盾或变暗结论。

两个公开 API:

* :func:`build_photometric_release` 消费同包 ``process_exposures`` 风格的
  测光批次,叠加历元、转移记录与灵敏度证明,产出带版本号、可 JSON 往返的
  测光发布件;
* :func:`assess_fade` 仅依据发布件(含经 JSON 往返、跨进程重建的发布件)
  回答“前帧可信实测 + 后帧清洁未检出”能否确认变暗到指定幅度。

发布件只携带派生观测证据与定标结构,不打包原始曝光行、坐标、匹配候选边表
或整份 photometry。
"""
from __future__ import annotations

from fractions import Fraction

from .checks import array, fields, ident, number
from .exposures import _decimal_str, _geometric_evidence, _validate

FORMAT = "star-crossmatch-photometric-release/1"

EXPOSURE_LIMIT = 32
REFERENCE_LIMIT = 256
DETECTION_LIMIT = 2000
TRANSFER_LIMIT = 128
LIMIT_ROW_LIMIT = 5000

# 帧级绝对零点状态。
ZP_PUBLISHED = "published"
ZP_INDETERMINATE = "indeterminate"
ZP_UNSOLVABLE = "unsolvable"
ZP_UNANCHORED = "unanchored"

# 分量级状态。
COMP_SOLVABLE = "solvable"
COMP_UNSOLVABLE = "unsolvable"

# assess_fade 结果状态。
FADE_CONFIRMED = "fade_confirmed"
FADE_INSUFFICIENT = "insufficient_evidence"
NOT_COMPARABLE = "not_comparable"
CALIBRATION_UNSOLVABLE = "calibration_unsolvable"
NO_MEASUREMENT = "no_trusted_measurement"
NO_SENSITIVITY_PROOF = "no_sensitivity_proof"
PROOF_BLOCKED = "sensitivity_proof_blocked"
INVALID_QUERY = "invalid_query"

# 虚拟“绝对零点”节点:把绝对基准统一成 value[曝光] - value[_ZERO] = 零点。
# 合法 ID 字符集为 [A-Za-z0-9_:-],哨兵含点号,不可能与曝光身份冲突。
_ZERO = "__absolute_zero__.sentinel"


# ---------------------------------------------------------------------------
# 带权(带偏移)并查集
# ---------------------------------------------------------------------------


class _OffsetDSU:
    """以 ``value[x] - value[parent[x]] = diff[x]`` 维护的带偏移并查集。

    不做路径压缩(本批规模极小),逻辑直白且全程精确有理:每次 find 沿父链
    累积偏移,返回 (根, value[x]-value[根])。
    """

    def __init__(self, nodes):
        self.parent = {x: x for x in nodes}
        self.diff = {x: Fraction(0) for x in nodes}

    def find(self, x):
        total = Fraction(0)
        while self.parent[x] != x:
            total += self.diff[x]
            x = self.parent[x]
        return x, total

    def union(self, a, b, delta):
        """加入约束 value[b] - value[a] = delta;返回 False 表示矛盾。"""
        ra, da = self.find(a)  # value[a] = value[ra] + da
        rb, db = self.find(b)  # value[b] = value[rb] + db
        if ra == rb:
            return db - da == delta
        # value[rb] - value[ra] = da + delta - db,把 rb 挂到 ra 下。
        self.parent[rb] = ra
        self.diff[rb] = da + delta - db
        return True

    def anchor(self, node, value, zero):
        """加入约束 value[node] = value,即 value[node] - value[zero] = value。"""
        return self.union(zero, node, value)

    def determinate(self, node, zero):
        root, _ = self.find(node)
        return root == self.find(zero)[0]

    def value_of(self, node, zero):
        """节点的确定绝对值(value[node]-value[zero]);未连通时为 None。"""
        root, offset = self.find(node)
        zero_root, zero_offset = self.find(zero)
        if root != zero_root:
            return None
        return offset - zero_offset

    def connected(self, a, b):
        return self.find(a)[0] == self.find(b)[0]


# ---------------------------------------------------------------------------
# 输入校验
# ---------------------------------------------------------------------------


def _validate_photometry_size(photometry):
    """新规模上限先于原曝光校验执行;字段与类型随后由原校验负责。"""
    if not isinstance(photometry, dict):
        raise ValueError("photometry 必须是对象")
    exposures = photometry.get("exposures")
    reference = photometry.get("reference")
    if not isinstance(exposures, list) or not 1 <= len(exposures) <= EXPOSURE_LIMIT:
        raise ValueError(f"exposures 数量必须为 1..{EXPOSURE_LIMIT}")
    if not isinstance(reference, list) or not 1 <= len(reference) <= REFERENCE_LIMIT:
        raise ValueError(f"reference 数量必须为 1..{REFERENCE_LIMIT}")
    total = 0
    seen_exp = set()
    for exp in exposures:
        if not isinstance(exp, dict) or not isinstance(exp.get("observations"), list):
            raise ValueError("exposure 必须含 observations 数组")
        total += len(exp["observations"])
        eid = exp.get("id")
        if eid is not None:
            if eid in seen_exp:
                raise ValueError("duplicate exposure")
            seen_exp.add(eid)
    if not 1 <= total <= DETECTION_LIMIT:
        raise ValueError(f"全批检测总数必须为 1..{DETECTION_LIMIT}")


def _validate_epochs(epochs, exposure_ids):
    if not isinstance(epochs, list) or not 1 <= len(epochs) <= EXPOSURE_LIMIT:
        raise ValueError(f"epochs 必须为 1..{EXPOSURE_LIMIT} 行的数组")
    rows = []
    seen = set()
    for row in epochs:
        fields(row, ("exposure", "mjd"))
        eid = ident(row["exposure"])
        if eid not in exposure_ids:
            raise ValueError(f"epochs 引用未知曝光: {eid}")
        if eid in seen:
            raise ValueError(f"epochs 对曝光 {eid} 重复登记")
        seen.add(eid)
        mjd = number(row["mjd"])
        if mjd < 0:
            raise ValueError("MJD 必须为非负十进制天数")
        rows.append((eid, mjd))
    if seen != exposure_ids:
        raise ValueError("每次曝光必须恰有一行 epochs")
    rows.sort(key=lambda x: x[0])
    return rows


def _validate_transfers(transfers, exposure_ids):
    if not isinstance(transfers, list) or not 0 <= len(transfers) <= TRANSFER_LIMIT:
        raise ValueError(f"transfers 必须为 0..{TRANSFER_LIMIT} 行的数组")
    rows = []
    seen_ids = set()
    for row in transfers:
        fields(row, ("id", "from_exposure", "to_exposure", "delta"))
        tid = ident(row["id"])
        if tid in seen_ids:
            raise ValueError(f"transfer id 批内重复: {tid}")
        seen_ids.add(tid)
        src = ident(row["from_exposure"])
        dst = ident(row["to_exposure"])
        if src not in exposure_ids or dst not in exposure_ids:
            raise ValueError(f"transfer {tid} 引用未知曝光")
        if src == dst:
            raise ValueError(f"transfer {tid} 两端必须是不同曝光")
        delta = number(row["delta"])  # value[to] - value[from] = delta
        rows.append({"id": tid, "from": src, "to": dst, "delta": delta})
    rows.sort(key=lambda x: x["id"])
    return rows


def _validate_limits(limits, exposure_ids, reference_ids):
    if not isinstance(limits, list) or not 0 <= len(limits) <= LIMIT_ROW_LIMIT:
        raise ValueError(f"limits 必须为 0..{LIMIT_ROW_LIMIT} 行的数组")
    rows = []
    seen = set()
    for row in limits:
        fields(row, ("exposure", "reference", "limiting_instrumental_magnitude"))
        eid = ident(row["exposure"])
        rid = ident(row["reference"])
        if eid not in exposure_ids:
            raise ValueError(f"limits 引用未知曝光: {eid}")
        if rid not in reference_ids:
            raise ValueError(f"limits 引用未知参考源: {rid}")
        key = (eid, rid)
        if key in seen:
            raise ValueError(f"limits 对 {eid}/{rid} 重复登记")
        seen.add(key)
        bound = number(row["limiting_instrumental_magnitude"])
        rows.append({"exposure": eid, "reference": rid, "limit": bound})
    rows.sort(key=lambda x: (x["exposure"], x["reference"]))
    return rows


def _validate_top(data):
    if not isinstance(data, dict):
        raise ValueError("输入必须是对象")
    fields(data, ("photometry", "epochs", "transfers", "limits"))


# ---------------------------------------------------------------------------
# 定标分量与假设
# ---------------------------------------------------------------------------


def _undirected_components(exposure_ids, transfers):
    """按全部转移记录(忽略方向)划分无向连通分量,无边曝光自成一组。"""
    adjacency = {eid: set() for eid in exposure_ids}
    for tr in transfers:
        adjacency[tr["from"]].add(tr["to"])
        adjacency[tr["to"]].add(tr["from"])
    seen = set()
    groups = []
    for start in exposure_ids:
        if start in seen:
            continue
        stack = [start]
        members = []
        seen.add(start)
        while stack:
            node = stack.pop()
            members.append(node)
            for nxt in adjacency[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        groups.append(sorted(members))
    # 分量身份按成员最小曝光 ID(等价于按成员元组)排序,与输入行序无关。
    groups.sort(key=lambda members: members)
    return groups


def _evaluate_hypothesis(members, edges, anchors, zero):
    """在给定转移子集与不可删除的绝对基准下评估一个假设。

    edges 为本假设保留的转移记录;anchors 为 {曝光: 绝对零点 Fraction}。
    返回 (可行, dsu);dsu 同时刻画假设内仍互通的子组与确定绝对值。
    """
    dsu = _OffsetDSU(set(members) | {zero})
    for tr in edges:
        if not dsu.union(tr["from"], tr["to"], tr["delta"]):
            return False, dsu
    for eid, value in anchors.items():
        if eid in members and not dsu.anchor(eid, value, zero):
            return False, dsu
    return True, dsu


def _component_hypotheses(members, transfers, anchors, zero):
    """枚举一个分量的全部可行假设(不删或恰好删一条)。

    返回按 (excluded is None 优先, 其次 excluded ID) 排序的假设列表,每项
    为 {"excluded": id|None, "edges": [...], "dsu": dsu}。无可行假设时为空。
    """
    local = [tr for tr in transfers if tr["from"] in members]
    hypotheses = []
    ok, dsu = _evaluate_hypothesis(members, local, anchors, zero)
    if ok:
        # 全部约束同时成立:只保留“不删记录”这一种,不再枚举单条排除。
        return [{"excluded": None, "edges": list(local), "dsu": dsu}]
    for tr in local:
        edges = [edge for edge in local if edge["id"] != tr["id"]]
        ok, dsu = _evaluate_hypothesis(members, edges, anchors, zero)
        if ok:
            hypotheses.append(
                {"excluded": tr["id"], "edges": edges, "dsu": dsu}
            )
    hypotheses.sort(key=lambda h: ("" if h["excluded"] is None else h["excluded"]))
    return hypotheses


# ---------------------------------------------------------------------------
# 发布件构建
# ---------------------------------------------------------------------------


def build_photometric_release(data):
    """构建可离线交接、可 JSON 往返的测光发布件。

    输入 ``data`` 恰含 photometry/epochs/transfers/limits 四项;校验失败或
    不支持的数据以 :class:`ValueError` 拒绝,不修改调用者数据,不交付部分
    结果。返回 :class:`PhotometricRelease`(dict 子类),可直接 ``json.dumps``。
    """
    _validate_top(data)
    photometry = data["photometry"]
    _validate_photometry_size(photometry)
    # 原曝光输入校验照常执行(字段、类型、ID、坐标、语法与原规模上限)。
    threshold, tolerance, reference, exposures = _validate(photometry)
    exposure_ids = {exp["id"] for exp in exposures}
    reference_ids = set(reference)
    epoch_rows = _validate_epochs(data["epochs"], exposure_ids)
    transfer_rows = _validate_transfers(data["transfers"], exposure_ids)
    limit_rows = _validate_limits(data["limits"], exposure_ids, reference_ids)

    # 几何证据:身份只看几何最优方案,锚点/零点沿用原锚点规则。
    evidence = _geometric_evidence(threshold, tolerance, reference, exposures)
    evidence_by_exp = {item["exposure"]: item for item in evidence}

    anchors = {}
    for eid in exposure_ids:
        zp = evidence_by_exp[eid]["zero_point"]
        if zp is not None:
            anchors[eid] = zp

    # 可信实测:参考源在所有几何最优方案中都配到同一条检测。
    measurements = []
    for eid in sorted(exposure_ids):
        forced = evidence_by_exp[eid]["forced"]
        for rid in sorted(forced):
            obs_id, instr = forced[rid]
            measurements.append(
                {
                    "exposure": eid,
                    "reference": rid,
                    "observation": obs_id,
                    "instrumental_magnitude": _decimal_str(instr),
                }
            )

    # 灵敏度证明:仅当参考源未出现在该帧任一检测的原始几何候选中时才清洁;
    # 被阻止的组合保留命中检测 ID,使拒绝结论可追溯。
    proofs = []
    for row in limit_rows:
        eid, rid = row["exposure"], row["reference"]
        blocking = sorted(evidence_by_exp[eid]["incidence"].get(rid, set()))
        proofs.append(
            {
                "exposure": eid,
                "reference": rid,
                "limiting_instrumental_magnitude": _decimal_str(row["limit"]),
                "usable": not blocking,
                "blocking_detections": blocking,
            }
        )

    # 几何身份派生记录:仅保留“有几何候选但未在所有几何最优方案中锁定同一
    # 检测”的(曝光,参考源)组合及命中检测 ID,使无可信实测的拒绝结论可追溯。
    # 不含坐标、角距或候选边表;已锁定为可信实测的组合不重复出现。
    geometric_identities = []
    forced_keys = {(m["exposure"], m["reference"]) for m in measurements}
    for eid in sorted(exposure_ids):
        for rid in sorted(evidence_by_exp[eid]["incidence"]):
            if (eid, rid) in forced_keys:
                continue
            geometric_identities.append(
                {
                    "exposure": eid,
                    "reference": rid,
                    "candidate_detections": sorted(
                        evidence_by_exp[eid]["incidence"][rid]
                    ),
                }
            )

    groups = _undirected_components(sorted(exposure_ids), transfer_rows)
    components_payload = []
    model_components = []
    exposure_status = {}
    exposure_value = {}
    zero = _ZERO
    for cindex, members in enumerate(groups):
        cid = f"component-{cindex + 1}"
        hypotheses = _component_hypotheses(members, transfer_rows, anchors, zero)
        if not hypotheses:
            status = COMP_UNSOLVABLE
            for eid in members:
                exposure_status[eid] = ZP_UNSOLVABLE
                exposure_value[eid] = None
        else:
            status = COMP_SOLVABLE
            values_per_exp = {eid: set() for eid in members}
            for h in hypotheses:
                for eid in members:
                    value = h["dsu"].value_of(eid, zero)
                    if value is not None:
                        values_per_exp[eid].add(value)
            for eid in members:
                values = values_per_exp[eid]
                if len(values) == 1:
                    exposure_status[eid] = ZP_PUBLISHED
                    exposure_value[eid] = _decimal_str(next(iter(values)))
                elif values:
                    exposure_status[eid] = ZP_INDETERMINATE
                    exposure_value[eid] = None
                else:
                    exposure_status[eid] = ZP_UNANCHORED
                    exposure_value[eid] = None
        hyp_payload = [
            {
                "hypothesis": f"{cid}#h{hindex + 1}",
                "excluded_transfer": h["excluded"],
            }
            for hindex, h in enumerate(hypotheses)
        ]
        components_payload.append(
            {
                "component": cid,
                "exposures": members,
                "status": status,
                "hypotheses": hyp_payload,
            }
        )
        model_components.append(
            {
                "component": cid,
                "exposures": list(members),
                "status": status,
                "hypotheses": hypotheses,
            }
        )

    exposures_payload = []
    for eid in sorted(exposure_ids):
        exposures_payload.append(
            {
                "exposure": eid,
                "mjd": _decimal_str(dict(epoch_rows)[eid]),
                "is_anchor": eid in anchors,
                "zero_point_status": exposure_status[eid],
                "zero_point": exposure_value[eid],
            }
        )

    transfers_payload = [
        {
            "id": tr["id"],
            "from_exposure": tr["from"],
            "to_exposure": tr["to"],
            "delta": _decimal_str(tr["delta"]),
            "component": next(
                c["component"]
                for c in model_components
                if tr["from"] in c["exposures"]
            ),
        }
        for tr in transfer_rows
    ]

    payload = {
        "format": FORMAT,
        "references": sorted(reference_ids),
        "exposures": exposures_payload,
        "transfers": transfers_payload,
        "components": components_payload,
        "measurements": measurements,
        "sensitivity_proofs": proofs,
        "geometric_identities": geometric_identities,
    }
    model = _build_model(payload)
    return PhotometricRelease(payload, model)


# ---------------------------------------------------------------------------
# 发布件装载与内存模型
# ---------------------------------------------------------------------------


def _build_model(payload):
    """从纯 JSON 结构重建精确有理内存模型(新进程无任何原始输入依赖)。"""
    if not isinstance(payload, dict) or payload.get("format") != FORMAT:
        raise ValueError("不支持的测光发布件格式")
    for key in (
        "references",
        "exposures",
        "transfers",
        "components",
        "measurements",
        "sensitivity_proofs",
        "geometric_identities",
    ):
        if key not in payload:
            raise ValueError(f"发布件缺少字段: {key}")

    references = list(payload["references"])
    exp_info = {}
    exp_order = []
    for row in payload["exposures"]:
        fields(
            row,
            ("exposure", "mjd", "is_anchor", "zero_point_status", "zero_point"),
        )
        eid = ident(row["exposure"])
        if eid in exp_info:
            raise ValueError("发布件曝光身份重复")
        if type(row["is_anchor"]) is not bool:
            raise ValueError("is_anchor 必须为布尔值")
        mjd = number(row["mjd"])
        if mjd < 0:
            raise ValueError("发布件 MJD 非法")
        zp = None if row["zero_point"] is None else number(row["zero_point"])
        exp_info[eid] = {
            "mjd": mjd,
            "is_anchor": row["is_anchor"],
            "zero_point_status": row["zero_point_status"],
            "zero_point": zp,
        }
        exp_order.append(eid)

    transfers = []
    for row in payload["transfers"]:
        fields(
            row,
            ("id", "from_exposure", "to_exposure", "delta", "component"),
        )
        transfers.append(
            {
                "id": ident(row["id"]),
                "from": ident(row["from_exposure"]),
                "to": ident(row["to_exposure"]),
                "delta": number(row["delta"]),
                "component": row["component"],
            }
        )

    components = []
    for crow in payload["components"]:
        fields(crow, ("component", "exposures", "status", "hypotheses"))
        cid = crow["component"]
        members = [ident(x) for x in crow["exposures"]]
        if any(eid not in exp_info for eid in members):
            raise ValueError(f"分量 {cid} 含未知曝光")
        local_edges = [
            tr for tr in transfers if tr["from"] in members and tr["to"] in members
        ]
        local_ids = {tr["id"] for tr in local_edges}
        anchors = {
            eid: exp_info[eid]["zero_point"]
            for eid in members
            if exp_info[eid]["is_anchor"]
        }
        zero = _ZERO
        hypotheses = []
        excluded_seen = set()
        for hrow in crow["hypotheses"]:
            fields(hrow, ("hypothesis", "excluded_transfer"))
            excluded = hrow["excluded_transfer"]
            if excluded is not None:
                excluded = ident(excluded)
                if excluded not in local_ids:
                    raise ValueError(f"假设排除了不属于 {cid} 的转移")
            if excluded in excluded_seen:
                raise ValueError(f"{cid} 假设身份重复")
            excluded_seen.add(excluded)
            edges = (
                list(local_edges)
                if excluded is None
                else [tr for tr in local_edges if tr["id"] != excluded]
            )
            ok, dsu = _evaluate_hypothesis(members, edges, anchors, zero)
            if not ok:
                raise ValueError(f"发布件 {cid} 含不可行假设")
            hypotheses.append(
                {
                    "hypothesis": hrow["hypothesis"],
                    "excluded": excluded,
                    "edges": edges,
                    "dsu": dsu,
                }
            )
        if crow["status"] == COMP_UNSOLVABLE and hypotheses:
            raise ValueError(f"{cid} 状态与假设矛盾")
        if crow["status"] == COMP_SOLVABLE and not hypotheses:
            raise ValueError(f"{cid} 状态与假设矛盾")
        components.append(
            {
                "component": cid,
                "exposures": members,
                "status": crow["status"],
                "hypotheses": hypotheses,
            }
        )

    component_of = {}
    for comp in components:
        for eid in comp["exposures"]:
            component_of[eid] = comp["component"]
    for tr in transfers:
        if component_of.get(tr["from"]) != tr["component"]:
            raise ValueError(f"转移 {tr['id']} 的分量归属不一致")

    measurements = {}
    for row in payload["measurements"]:
        fields(
            row,
            ("exposure", "reference", "observation", "instrumental_magnitude"),
        )
        eid, rid = ident(row["exposure"]), ident(row["reference"])
        if eid not in exp_info or rid not in references:
            raise ValueError("发布件实测含未知身份")
        key = (eid, rid)
        if key in measurements:
            raise ValueError("发布件实测身份重复")
        measurements[key] = {
            "observation": ident(row["observation"]),
            "instrumental": number(row["instrumental_magnitude"]),
        }

    proofs = {}
    for row in payload["sensitivity_proofs"]:
        fields(
            row,
            (
                "exposure",
                "reference",
                "limiting_instrumental_magnitude",
                "usable",
                "blocking_detections",
            ),
        )
        eid, rid = ident(row["exposure"]), ident(row["reference"])
        if eid not in exp_info or rid not in references:
            raise ValueError("发布件灵敏度证明含未知身份")
        if type(row["usable"]) is not bool:
            raise ValueError("usable 必须为布尔值")
        key = (eid, rid)
        if key in proofs:
            raise ValueError("发布件灵敏度证明身份重复")
        blocking = [ident(x) for x in row["blocking_detections"]]
        if row["usable"] and blocking:
            raise ValueError("清洁证明不得附带阻止检测")
        if not row["usable"] and not blocking:
            raise ValueError("被阻止证明必须给出检测 ID")
        proofs[key] = {
            "limit": number(row["limiting_instrumental_magnitude"]),
            "usable": row["usable"],
            "blocking_detections": sorted(blocking),
        }

    geometric = {}
    for row in payload["geometric_identities"]:
        fields(row, ("exposure", "reference", "candidate_detections"))
        eid, rid = ident(row["exposure"]), ident(row["reference"])
        if eid not in exp_info or rid not in references:
            raise ValueError("发布件几何身份记录含未知身份")
        key = (eid, rid)
        if key in geometric:
            raise ValueError("发布件几何身份记录重复")
        if key in measurements:
            raise ValueError("已锁定可信实测的组合不得再登记几何身份记录")
        candidate_ids = sorted(ident(x) for x in row["candidate_detections"])
        if not candidate_ids:
            raise ValueError("几何身份记录必须至少含一条候选检测")
        geometric[key] = candidate_ids

    return {
        "references": references,
        "exposures": exp_info,
        "exposure_order": exp_order,
        "components": components,
        "component_of": component_of,
        "measurements": measurements,
        "proofs": proofs,
        "geometric": geometric,
    }


class PhotometricRelease(dict):
    """测光发布件;dict 本身即版本化 JSON 载荷,另持有精确有理内存模型。

    经 ``json.dumps``/``json.loads`` 往返后,把纯 dict 传给
    :func:`assess_fade` 即可,查询结论与原进程完全一致。
    """

    def __init__(self, payload, model=None):
        super().__init__(payload)
        self._model = model if model is not None else _build_model(payload)

    @classmethod
    def from_dict(cls, payload):
        """从(可能经 JSON 往返的)纯 dict 重建发布件并严格校验。"""
        if isinstance(payload, PhotometricRelease):
            return cls(dict(payload), payload._model)
        return cls(payload, _build_model(payload))

    def to_json(self, **kwargs):
        """序列化为 JSON 文本(默认 sort_keys、拒绝非有限值)。"""
        import json

        options = {"ensure_ascii": False, "sort_keys": True, "allow_nan": False}
        options.update(kwargs)
        return json.dumps(dict(self), **options)


# ---------------------------------------------------------------------------
# 变暗判据
# ---------------------------------------------------------------------------


def _validate_query(query):
    if not isinstance(query, dict):
        raise ValueError("query 必须是对象")
    fields(
        query,
        ("reference", "bright_exposure", "later_exposure", "minimum_fade"),
    )
    rid = ident(query["reference"])
    bright = ident(query["bright_exposure"])
    later = ident(query["later_exposure"])
    minimum = number(query["minimum_fade"])
    if minimum <= 0:
        raise ValueError("minimum_fade 必须为正十进制星等")
    return rid, bright, later, minimum


def _evidence_ref(component_id, excluded):
    return {"component": component_id, "excluded_transfer": excluded}


def assess_fade(release, query):
    """依据发布件复核一次变暗查询,返回可 JSON 序列化的结论 dict。

    结论状态见模块级常量;物理上的定标冲突与证据不足是正常业务结果,
    只有未知身份、坏数据、逆时间查询或不支持的发布件格式才抛
    :class:`ValueError`。
    """
    if isinstance(release, PhotometricRelease):
        model = release._model
    elif isinstance(release, dict):
        model = _build_model(release)
    else:
        raise ValueError("需要测光发布件 dict 或 PhotometricRelease")

    rid, bright, later, minimum = _validate_query(query)
    base = {
        "reference": rid,
        "bright_exposure": bright,
        "later_exposure": later,
        "minimum_fade": _decimal_str(minimum),
    }

    def result(status, **extra):
        out = dict(base)
        out["status"] = status
        out.update(extra)
        return out

    if rid not in model["references"]:
        raise ValueError(f"查询引用未知参考源: {rid}")
    if bright not in model["exposures"] or later not in model["exposures"]:
        raise ValueError("查询引用未知曝光")
    mjd_bright = model["exposures"][bright]["mjd"]
    mjd_later = model["exposures"][later]["mjd"]
    if not mjd_bright < mjd_later:
        raise ValueError("bright_exposure 的 MJD 必须严格小于 later_exposure")

    # 1. 早于定标判断的三类观测证据拒绝,顺序固定、状态可区分。
    measurement = model["measurements"].get((bright, rid))
    if measurement is None:
        return result(
            NO_MEASUREMENT,
            detections=list(model["geometric"].get((bright, rid), [])),
            detail="前帧没有覆盖全部几何最优方案的可信实测",
        )
    proof = model["proofs"].get((later, rid))
    if proof is None:
        return result(
            NO_SENSITIVITY_PROOF,
            detections=[],
            detail="后帧没有该源的灵敏度证明",
        )
    if not proof["usable"]:
        return result(
            PROOF_BLOCKED,
            detections=list(proof["blocking_detections"]),
            detail="后帧存在该源的原始几何候选,未检出证明不可用",
        )

    # 2. 定标可行性与联合假设。
    comp_bright = model["component_of"][bright]
    comp_later = model["component_of"][later]
    comps = {c["component"]: c for c in model["components"]}
    involved = []
    for cid in (comp_bright, comp_later):
        if cid not in involved:
            involved.append(cid)
    for cid in involved:
        if comps[cid]["status"] == COMP_UNSOLVABLE:
            return result(
                CALIBRATION_UNSOLVABLE,
                unsolvable_components=list(involved),
                detail="所用分量没有可行假设",
            )

    # 同一原分量的两帧必须在同一个排除假设下作差;不同分量的假设任意组合。
    # 联合假设按分量 ID、再按假设在发布件中的固定次序枚举,首个即字典序最小,
    # 作为“并列”文档固定规则下的证据假设。
    joint_options = []
    if comp_bright == comp_later:
        comp = comps[comp_bright]
        for hindex, hyp in enumerate(comp["hypotheses"]):
            joint_options.append(
                {
                    "indices": {comp_bright: hindex},
                    "evidence": [
                        _evidence_ref(comp_bright, hyp["excluded"])
                    ],
                }
            )
    else:
        cb, cl = comps[comp_bright], comps[comp_later]
        for hi, hb in enumerate(cb["hypotheses"]):
            for hj, hl in enumerate(cl["hypotheses"]):
                joint_options.append(
                    {
                        "indices": {comp_bright: hi, comp_later: hj},
                        "evidence": [
                            _evidence_ref(comp_bright, hb["excluded"]),
                            _evidence_ref(comp_later, hl["excluded"]),
                        ],
                    }
                )

    zero = _ZERO
    instr_meas = measurement["instrumental"]
    limit = proof["limit"]
    bounds = []
    first_incomparable = None
    for joint in joint_options:
        # 取该联合假设下两帧的零点(或确定绝对值)与互通关系。
        if comp_bright == comp_later:
            comp = comps[comp_bright]
            hyp = comp["hypotheses"][joint["indices"][comp_bright]]
            dsu = hyp["dsu"]
            connected = dsu.connected(bright, later)
            zb = dsu.value_of(bright, zero)
            zl = dsu.value_of(later, zero)
        else:
            hb = comps[comp_bright]["hypotheses"][
                joint["indices"][comp_bright]
            ]
            hl = comps[comp_later]["hypotheses"][
                joint["indices"][comp_later]
            ]
            connected = False
            zb = hb["dsu"].value_of(bright, zero)
            zl = hl["dsu"].value_of(later, zero)
        if not connected and (zb is None or zl is None):
            # 零点差不受约束:整份查询不可比较,返回字典序最小的对应假设。
            if first_incomparable is None:
                first_incomparable = joint
            continue
        if connected:
            # 同一互通子组内零点差精确确定(与是否有绝对基准无关)。
            _, off_b = dsu.find(bright)
            _, off_l = dsu.find(later)
            delta_zp = off_l - off_b
        else:
            delta_zp = zl - zb
        bound = limit - instr_meas + delta_zp
        bounds.append((bound, joint))

    if first_incomparable is not None:
        return result(
            NOT_COMPARABLE,
            hypothesis_evidence=first_incomparable["evidence"],
            detail="存在两帧不互通且至少一方无绝对基准的联合假设",
        )

    # 3. 全部联合假设都可比:取可保证的变暗下界。
    bounds.sort(key=lambda pair: pair[0])
    minimum_bound, winning_joint = bounds[0]
    confirmed = minimum_bound >= minimum
    return result(
        FADE_CONFIRMED if confirmed else FADE_INSUFFICIENT,
        guaranteed_fade_lower_bound=_decimal_str(minimum_bound),
        hypothesis_evidence=winning_joint["evidence"],
        detail=(
            "下界达到 minimum_fade,确认变暗"
            if confirmed
            else "下界低于 minimum_fade,证据不足,不表示源未变暗"
        ),
    )
