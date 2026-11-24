"""离线跟拍工作簿。

本模块只依赖标准库。工作簿把测光批次、当前待执行指令和已执行记录放在同一个
JSON 文件中；每次成功修改都通过临时文件 + ``os.replace`` 原子落盘，因此调用方
可以关闭后用 :meth:`FollowUpWorkbook.open` 重开。

所有统计判定均使用 :class:`fractions.Fraction`，不使用浮点结果决定候选资格。
"""

from __future__ import annotations

import copy
import json
import os
import tempfile
from math import comb
from pathlib import Path
from fractions import Fraction

from .checks import ident, integer, number
from .io import load as _load_json

FORMAT = "star-crossmatch-followup/1"
DEFAULT_WORKBOOK_NAME = "followup_workbook.json"


class FollowUpError(ValueError):
    """跟拍工作簿输入或状态冲突；同时是 ``ValueError``。"""


class FollowUpConflict(FollowUpError):
    """指令令牌或实拍记录身份冲突。"""


class AttributeDict(dict):
    """同时支持字典键和属性访问的结果对象。"""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


class CorrectionResult(AttributeDict):
    """整批更正结果；可通过 ``revoked``/``issued`` 键或属性访问。"""


# ---------------------------------------------------------------------------
# 通用小工具
# ---------------------------------------------------------------------------


class IndexedList(list):
    """既可按下标/迭代访问，也可用业务 ID 字符串查找的列表。"""

    def __init__(self, items, *key_fields):
        super().__init__(items)
        self._maps = []
        for field in key_fields:
            mapping = {}
            for item in items:
                if isinstance(item, dict) and field in item:
                    mapping[item[field]] = item
            self._maps.append((field, mapping))

    def __getitem__(self, index):
        if isinstance(index, str):
            for _field, mapping in self._maps:
                if index in mapping:
                    return mapping[index]
            raise KeyError(index)
        return list.__getitem__(self, index)

    def __getattr__(self, name):
        # 仅在普通属性查找失败后调用，因此不会覆盖 list/dict 风格方法。
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def get(self, key, default=None):
        if isinstance(key, str):
            for _field, mapping in self._maps:
                if key in mapping:
                    return mapping[key]
            return default
        return list.__getitem__(self, key)

    def __contains__(self, key):
        if isinstance(key, str):
            return any(key in mapping for _field, mapping in self._maps)
        return list.__contains__(self, key)

    def keys(self):
        return self._maps[0][1].keys()

    def values(self):
        return self._maps[0][1].values()

    def items(self):
        return self._maps[0][1].items()


def _strict_bool(value, what="布尔标记"):
    if type(value) is not bool:
        raise FollowUpError(f"{what}必须是 JSON 布尔值")
    return value


def _fraction(value, what="十进制字符串", *, nonnegative=False, positive=False):
    if isinstance(value, Fraction):
        result = value
    elif type(value) is int:
        result = Fraction(value)
    elif isinstance(value, str):
        try:
            result = number(value)
        except ValueError as exc:
            raise FollowUpError(f"{what}不是仓库语法允许的十进制字符串") from exc
    else:
        raise FollowUpError(f"{what}必须是十进制字符串")
    if nonnegative and result < 0:
        raise FollowUpError(f"{what}必须非负")
    if positive and result <= 0:
        raise FollowUpError(f"{what}必须为正")
    return result


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _pick(row, names, what):
    present = [name for name in names if name in row]
    if not present:
        raise FollowUpError(f"{what}缺少字段，候选名称之一：{', '.join(names)}")
    values = [row[name] for name in present]
    if any(value != values[0] for value in values[1:]):
        raise FollowUpError(f"{what}的别名字段给出了冲突值")
    return values[0]


def _object(row, allowed, what):
    if not isinstance(row, dict):
        raise FollowUpError(f"{what}必须是对象")
    extra = set(row) - set(allowed)
    if extra:
        raise FollowUpError(f"{what}含未定义字段：{sorted(extra)}")
    return row


def _array(value, low, high, what):
    if not isinstance(value, list) or not low <= len(value) <= high:
        raise FollowUpError(f"{what}必须是 {low}..{high} 项的数组")
    return value


# ---------------------------------------------------------------------------
# 输入解析与持久化状态
# ---------------------------------------------------------------------------

_SOURCE_ID_KEYS = ("id", "source_id", "source")
_BASELINE_KEYS = ("baseline_magnitude", "baseline", "baseline_mag", "reference_magnitude")
_PRIORITY_KEYS = ("priority", "followup_priority", "follow_up_priority", "order")

_MEASUREMENT_ID_KEYS = ("id", "measurement_id", "measurement")
_MEASUREMENT_SOURCE_KEYS = ("source_id", "source")
_MJD_KEYS = ("mjd", "mjd_days", "time")
_MAGNITUDE_KEYS = ("magnitude", "mag")
_SIGMA_KEYS = ("sigma", "magnitude_sigma", "sigma_magnitude", "uncertainty")
_USABLE_KEYS = ("usable", "is_usable", "useful")

_TOP_SOURCES = ("sources",)
_TOP_MEASUREMENTS = ("measurements", "photometry", "observations")
_TOP_MIN_TIMES = (
    "min_valid_times",
    "minimum_valid_times",
    "min_times",
    "min_valid_count",
    "min_valid_time_points",
    "min_valid_epochs",
    "min_epochs",
)
_TOP_FDR = (
    "fdr_threshold",
    "false_discovery_rate_threshold",
    "q_threshold",
    "fdr_q",
    "fdr_alpha",
    "fdr",
    "threshold",
    "alpha",
)
_TOP_CAPACITY = (
    "capacity",
    "followup_capacity",
    "follow_up_capacity",
    "slot_count",
    "slots",
    "quota",
    "max_followups",
)

_TOP_FIELDS = set(
    _TOP_SOURCES
    + _TOP_MEASUREMENTS
    + _TOP_MIN_TIMES
    + _TOP_FDR
    + _TOP_CAPACITY
)
_SOURCE_FIELDS = set(_SOURCE_ID_KEYS + _BASELINE_KEYS + _PRIORITY_KEYS)
_MEASUREMENT_FIELDS = set(
    _MEASUREMENT_ID_KEYS
    + _MEASUREMENT_SOURCE_KEYS
    + _MJD_KEYS
    + _MAGNITUDE_KEYS
    + _SIGMA_KEYS
    + _USABLE_KEYS
)


def _merge_batch(data, options):
    batch = {}
    if data is not None:
        if not isinstance(data, dict):
            raise FollowUpError("初始化数据必须是批次对象；源列表可通过 sources= 传入")
        batch.update(data)
    for key, value in options.items():
        if value is not None:
            batch[key] = value
    extra = set(batch) - _TOP_FIELDS
    if extra:
        raise FollowUpError(f"初始化批次含未定义字段：{sorted(extra)}")
    return batch


def _build_state(batch: dict) -> dict:
    if not isinstance(batch, dict):
        raise FollowUpError("初始化批次必须是对象")

    source_rows = _array(
        _pick(batch, _TOP_SOURCES, "初始化批次"), 1, 1000, "sources"
    )
    measurement_rows = _array(
        _pick(batch, _TOP_MEASUREMENTS, "初始化批次"),
        1,
        10000,
        "measurements",
    )
    min_valid_times = integer(
        _pick(batch, _TOP_MIN_TIMES, "初始化批次"), 3, 128
    )
    threshold = _fraction(
        _pick(batch, _TOP_FDR, "初始化批次"), "误发现率阈值"
    )
    if not 0 < threshold <= 1:
        raise FollowUpError("误发现率阈值必须位于 (0,1]")

    # 第一遍：源身份、基准星等和优先序。
    sources = []
    source_ids = set()
    priorities = []
    for index, row in enumerate(source_rows):
        _object(row, _SOURCE_FIELDS, f"sources[{index}]")
        sid = ident(_pick(row, _SOURCE_ID_KEYS, f"sources[{index}]"))
        if sid in source_ids:
            raise FollowUpError(f"重复源 ID：{sid}")
        source_ids.add(sid)

        baseline = _fraction(
            _pick(row, _BASELINE_KEYS, f"sources[{index}]"), "基准星等"
        )
        priority = integer(
            _pick(row, _PRIORITY_KEYS, f"sources[{index}]"), 1, len(source_rows)
        )
        priorities.append(priority)
        sources.append(
            {
                "id": sid,
                "baseline": _fraction_text(baseline),
                "priority": priority,
            }
        )

    source_count = len(sources)
    if sorted(priorities) != list(range(1, source_count + 1)):
        raise FollowUpError("跟拍优先序必须是 1..源数 的全批唯一排列")
    capacity = integer(
        _pick(batch, _TOP_CAPACITY, "初始化批次"), 1, source_count
    )

    # 第二遍：测量身份、已登记源、范围和每源不同时刻上限。
    measurements = []
    measurement_ids = set()
    times_by_source = {sid: set() for sid in source_ids}
    for index, row in enumerate(measurement_rows):
        _object(row, _MEASUREMENT_FIELDS, f"measurements[{index}]")
        mid = ident(_pick(row, _MEASUREMENT_ID_KEYS, f"measurements[{index}]"))
        if mid in measurement_ids:
            raise FollowUpError(f"重复测量 ID：{mid}")
        measurement_ids.add(mid)

        sid = ident(
            _pick(row, _MEASUREMENT_SOURCE_KEYS, f"measurements[{index}]")
        )
        if sid not in source_ids:
            raise FollowUpError(f"测量 {mid} 指向未登记源：{sid}")

        mjd = _fraction(
            _pick(row, _MJD_KEYS, f"measurements[{index}]"),
            "MJD",
            nonnegative=True,
        )
        magnitude = _fraction(
            _pick(row, _MAGNITUDE_KEYS, f"measurements[{index}]"), "星等"
        )
        sigma = _fraction(
            _pick(row, _SIGMA_KEYS, f"measurements[{index}]"),
            "星等标准差",
            positive=True,
        )
        usable = _strict_bool(
            _pick(row, _USABLE_KEYS, f"measurements[{index}]"), "usable"
        )

        times_by_source[sid].add(mjd)
        measurements.append(
            {
                "id": mid,
                "source_id": sid,
                "mjd": _fraction_text(mjd),
                "magnitude": _fraction_text(magnitude),
                "sigma": _fraction_text(sigma),
                "usable": usable,
            }
        )

    for sid, times in times_by_source.items():
        if len(times) > 128:
            raise FollowUpError(f"源 {sid} 有超过 128 个不同输入时刻")

    state = {
        "format": FORMAT,
        "config": {
            "min_valid_times": min_valid_times,
            "fdr_threshold": _fraction_text(threshold),
            "capacity": capacity,
        },
        "sources": sources,
        "measurements": measurements,
        "next_token": 1,
        "pending": [],
        "executions": [],
    }
    assessment = _assess_state(state)
    _reconcile_plan(state, assessment)
    return state


def _resolve_workbook_path(path, create=False):
    """把路径解析为 JSON 工作簿文件。

    若路径已经是目录，则使用目录下的固定文件名；否则路径本身就是文件，
    有没有 ``.json`` 后缀均可。
    """
    path = Path(path)
    if path.exists() and path.is_dir():
        return path / DEFAULT_WORKBOOK_NAME
    return path


def _fsync_directory(path: Path):
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    fd = os.open(path, flags)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _write_state(path: Path, state: dict, overwrite: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(
            state,
            ensure_ascii=False,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")

    temporary_name = None
    try:
        fd, temporary_name = tempfile.mkstemp(
            prefix=".followup-workbook-", suffix=".tmp", dir=path.parent
        )
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())

        if overwrite:
            os.replace(temporary_name, path)
            temporary_name = None
        else:
            # 硬链接到目标可以原子地保证“不覆盖既有工作簿”。
            os.link(temporary_name, path)
        _fsync_directory(path.parent)
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def _read_state(path: Path) -> dict:
    state = _load_json(path)
    if not isinstance(state, dict) or state.get("format") != FORMAT:
        raise FollowUpError("不是受支持的跟拍工作簿文件")
    for key in (
        "config",
        "sources",
        "measurements",
        "next_token",
        "pending",
        "executions",
    ):
        if key not in state:
            raise FollowUpError(f"工作簿文件缺少状态字段：{key}")
    return state


# ---------------------------------------------------------------------------
# 精确统计
# ---------------------------------------------------------------------------


def sign_test_probability(n: int, k: int) -> Fraction:
    """n 次公平硬币中至少 k 次正面的精确概率。"""
    if not 0 <= k <= n:
        raise FollowUpError("偏亮时刻数 k 必须位于 0..n")
    numerator = sum(comb(n, j) for j in range(k, n + 1))
    return Fraction(numerator, 2**n)


def _merged_magnitude(records):
    weight_sum = Fraction(0)
    weighted_sum = Fraction(0)
    for record in records:
        sigma = Fraction(record["sigma"])
        weight = 1 / (sigma * sigma)
        weight_sum += weight
        weighted_sum += Fraction(record["magnitude"]) * weight
    return weighted_sum / weight_sum


def _assess_state(state: dict) -> dict:
    min_valid_times = state["config"]["min_valid_times"]
    threshold = Fraction(state["config"]["fdr_threshold"])

    measurements_by_source = {source["id"]: [] for source in state["sources"]}
    for measurement in state["measurements"]:
        if measurement["usable"]:
            measurements_by_source[measurement["source_id"]].append(measurement)

    results = {}
    eligible = []
    for source in state["sources"]:
        sid = source["id"]
        baseline = Fraction(source["baseline"])
        groups = {}
        for measurement in measurements_by_source[sid]:
            mjd = Fraction(measurement["mjd"])
            groups.setdefault(mjd, []).append(measurement)

        k = 0
        non_equal = 0
        for records in groups.values():
            merged = _merged_magnitude(records)
            if merged < baseline:
                k += 1
                non_equal += 1
            elif merged > baseline:
                non_equal += 1
            # 恰好等于基准的时刻不计入 n、k。

        result = {
            "source_id": sid,
            "priority": source["priority"],
            "n": non_equal,
            "k": k,
            "p": None,
            "q": None,
            "candidate": False,
            "status": "insufficient_evidence",
        }
        if non_equal >= min_valid_times:
            result["p"] = sign_test_probability(non_equal, k)
            result["status"] = "not_candidate"
            eligible.append(result)
        results[sid] = result

    # Benjamini-Hochberg：先按 (p, 源 ID) 确定排名；反向最小值保证同 p 同 q。
    eligible.sort(key=lambda item: (item["p"], item["source_id"]))
    total = len(eligible)
    running = Fraction(1)
    for reverse_index in range(total - 1, -1, -1):
        rank = reverse_index + 1
        item = eligible[reverse_index]
        raw_q = total * item["p"] / rank
        running = min(running, raw_q)
        item["q"] = min(Fraction(1), running)

    for item in eligible:
        if item["q"] <= threshold:
            item["candidate"] = True
            item["status"] = "candidate"

    return results


def _public_candidate(row: dict) -> dict:
    p = row["p"]
    q = row["q"]
    candidate = bool(row["candidate"])
    status = row["status"]
    return AttributeDict({
        "source_id": row["source_id"],
        "source": row["source_id"],
        "id": row["source_id"],  # 便捷别名
        "priority": row["priority"],
        "n": row["n"],
        "k": row["k"],
        "p": None if p is None else _fraction_text(p),
        "q": None if q is None else _fraction_text(q),
        "candidate": candidate,
        "is_candidate": candidate,  # 语义别名
        "insufficient_evidence": status == "insufficient_evidence",
        "status": status,
    })


def _public_instruction(instruction: dict) -> dict:
    result = dict(instruction)
    result["instruction_id"] = instruction["token"]
    result["source"] = instruction["source_id"]  # 与仓库既有 *source ID 输出习惯兼容
    # 把常用别名字段放在前面，便于人读。
    ordered = AttributeDict()
    for key in (
        "token",
        "instruction_id",
        "source_id",
        "source",
        "n",
        "k",
        "p",
        "q",
    ):
        if key in result:
            ordered[key] = result[key]
    for key, value in result.items():
        if key not in ordered:
            ordered[key] = value
    return ordered


def _public_execution(execution: dict) -> dict:
    record_id = execution["record_id"]
    return AttributeDict({
        "token": execution["token"],
        "source_id": execution["instruction"]["source_id"],
        "source": execution["instruction"]["source_id"],
        "record_id": record_id,
        "id": record_id,
        "shot_id": record_id,  # 实拍记录 ID 的语义别名
        "capture_id": record_id,
        "instruction": _public_instruction(execution["instruction"]),
    })


# ---------------------------------------------------------------------------
# 计划调和
# ---------------------------------------------------------------------------


def _reconcile_plan(state: dict, assessment: dict):
    capacity = state["config"]["capacity"]
    executed_sources = {
        execution["instruction"]["source_id"]
        for execution in state["executions"]
    }
    remaining_slots = capacity - len(state["executions"])
    if remaining_slots < 0:
        raise FollowUpError("已执行记录数超过容量")

    priority_by_source = {
        source["id"]: source["priority"] for source in state["sources"]
    }
    desired = [
        sid
        for sid, row in assessment.items()
        if row["candidate"] and sid not in executed_sources
    ]
    desired.sort(key=lambda sid: priority_by_source[sid])
    desired = desired[:remaining_slots]
    desired_set = set(desired)

    existing_by_source = {
        instruction["source_id"]: instruction
        for instruction in state["pending"]
    }

    revoked = []
    retained_pending = []
    for instruction in state["pending"]:
        if instruction["source_id"] in desired_set:
            retained_pending.append(instruction)
        else:
            revoked.append(instruction)

    issued = []
    next_pending = []
    for sid in desired:
        if sid in existing_by_source:
            instruction = existing_by_source[sid]
        else:
            row = assessment[sid]
            token_number = state["next_token"]
            state["next_token"] = token_number + 1
            instruction = {
                "token": f"fu-{token_number:012d}",
                "source_id": sid,
                "n": row["n"],
                "k": row["k"],
                "p": _fraction_text(row["p"]),
                "q": _fraction_text(row["q"]),
            }
            issued.append(instruction)
        next_pending.append(instruction)

    state["pending"] = next_pending
    return {"revoked": revoked, "issued": issued}


# ---------------------------------------------------------------------------
# 更正输入
# ---------------------------------------------------------------------------


def _normalize_corrections(corrections, kwargs):
    raw_items = []
    if corrections is not None:
        if isinstance(corrections, dict) and any(
            key in corrections for key in _USABLE_KEYS
        ) and any(key in corrections for key in _MEASUREMENT_ID_KEYS):
            # 单个 {measurement_id, usable} 对象。
            raw_items.append(
                (
                    _pick(corrections, _MEASUREMENT_ID_KEYS, "质量更正"),
                    _pick(corrections, _USABLE_KEYS, "质量更正"),
                )
            )
        elif isinstance(corrections, dict) and any(
            key in corrections
            for key in (
                "corrections",
                "updates",
                "changes",
                "flags",
                "quality_updates",
            )
        ):
            for key in (
                "corrections",
                "updates",
                "changes",
                "flags",
                "quality_updates",
            ):
                if key in corrections:
                    return _normalize_corrections(corrections[key], kwargs)
        elif isinstance(corrections, dict):
            raw_items.extend(corrections.items())
        elif isinstance(corrections, (list, tuple)):
            for item in corrections:
                if isinstance(item, dict):
                    mid = _pick(item, _MEASUREMENT_ID_KEYS, "质量更正")
                    usable = _pick(item, _USABLE_KEYS, "质量更正")
                    raw_items.append((mid, usable))
                elif isinstance(item, (list, tuple)) and len(item) == 2:
                    raw_items.append((item[0], item[1]))
                else:
                    raise FollowUpError(
                        "质量更正项必须是 {id, usable} 或 (id, usable)"
                    )
        else:
            raise FollowUpError("质量更正必须是映射或 {id, usable} 列表")

    # 允许 update_quality_flags(m1=False, m2=True) 形式。
    raw_items.extend(kwargs.items())

    seen = set()
    normalized = []
    for mid, value in raw_items:
        mid = ident(mid)
        if mid in seen:
            raise FollowUpError(f"同一批质量更正重复指向测量：{mid}")
        seen.add(mid)

        if isinstance(value, dict):
            usable = _strict_bool(
                _pick(value, _USABLE_KEYS, f"测量 {mid} 的更正值"), "usable"
            )
        else:
            usable = _strict_bool(value, f"测量 {mid} 的 usable")
        normalized.append((mid, usable))

    return normalized


# ---------------------------------------------------------------------------
# 工作簿
# ---------------------------------------------------------------------------


class FollowUpWorkbook:
    """一个目录/文件对应一个跟拍窗口的本地持久化工作簿。"""

    def __init__(self, path=None, data=None, **options):
        create = data is not None or bool(options)
        raw_path = Path(path) if path is not None else None
        if raw_path is None:
            raise FollowUpError("必须给出工作簿路径")

        if create:
            batch = _merge_batch(data, options)
            # 先完整解析，再创建目录/文件，保证无效初始化不留下工作簿状态。
            state = _build_state(batch)
            self.path = _resolve_workbook_path(raw_path, True)
            _write_state(self.path, state, overwrite=False)
            self._state = state
        else:
            self.path = _resolve_workbook_path(raw_path, False)
            self._state = _read_state(self.path)

        self._closed = False

    # ---- 打开/初始化 -----------------------------------------------------

    @classmethod
    def create(cls, path, data=None, measurements=None, **options):
        """创建并初始化工作簿。

        可传完整批次字典 ``create(path, batch)``，也可传
        ``create(path, sources, measurements, min_valid_times=..., ...)``。
        """
        if isinstance(data, list):
            options["sources"] = data
            if measurements is not None:
                options["measurements"] = measurements
            data = None
        elif measurements is not None:
            options["measurements"] = measurements
        return cls(path, data, **options)

    @classmethod
    def initialize(cls, path, data=None, measurements=None, **options):
        return cls.create(path, data, measurements, **options)

    @classmethod
    def init(cls, path, data=None, measurements=None, **options):
        return cls.create(path, data, measurements, **options)

    @classmethod
    def open(cls, path):
        return cls(path)

    @classmethod
    def load(cls, path):
        return cls.open(path)

    def close(self):
        # 所有成功修改均已即时 fsync + atomic replace；保留该方法用于上下文管理。
        self._closed = True

    def save(self):
        """显式把当前状态再原子写盘一次；普通修改无需手动调用。"""
        _write_state(self.path, self._state, overwrite=True)

    def reload(self):
        """从磁盘重新读取工作簿（不支持并发写，仅用于显式刷新）。"""
        self._state = _read_state(self.path)
        return self

    refresh = reload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    # ---- 事务 ------------------------------------------------------------

    def _commit(self, operation):
        new_state = copy.deepcopy(self._state)
        raw_result = operation(new_state)
        # 完整写入成功后才替换当前内存状态。
        _write_state(self.path, new_state, overwrite=True)
        self._state = new_state
        return copy.deepcopy(raw_result)

    # ---- 配置 ------------------------------------------------------------

    @property
    def capacity(self):
        return int(self._state["config"]["capacity"])

    @property
    def execution_count(self):
        return len(self._state["executions"])

    @property
    def used_capacity(self):
        return len(self._state["executions"])

    @property
    def remaining_capacity(self):
        return self.capacity - self.execution_count

    def capacity_status(self):
        return {
            "capacity": self.capacity,
            "executed": self.execution_count,
            "used_capacity": self.execution_count,
            "remaining_capacity": self.remaining_capacity,
            "pending_count": len(self._state["pending"]),
        }

    # ---- 当前统计与候选 --------------------------------------------------

    def _public_candidate_rows(self):
        assessment = _assess_state(self._state)
        rows = [
            _public_candidate(assessment[source["id"]])
            for source in sorted(
                self._state["sources"], key=lambda source: source["priority"]
            )
        ]
        return IndexedList(rows, "source_id", "source", "id")

    def candidate_results(self):
        """全部源的当前 n/k/p/q 与 candidate/status 标记，按优先序返回。"""
        return self._public_candidate_rows()

    def candidate_map(self):
        return {row["source_id"]: row for row in self.candidate_results()}

    def candidates_by_id(self):
        return self.candidate_map()

    def source_results(self):
        return self.candidate_results()

    def statistics(self):
        return self.candidate_results()

    def stats(self):
        return self.candidate_results()

    def candidates(self, current_only=False):
        rows = self.candidate_results()
        if current_only:
            rows = [row for row in rows if row["candidate"]]
            return IndexedList(rows, "source_id", "source", "id")
        return rows

    def all_candidates(self):
        return self.candidate_results()

    def current_candidates(self):
        return self.candidates(current_only=True)

    list_current_candidates = current_candidates

    def get_source(self, source_id):
        source_id = ident(source_id)
        return self.candidate_results()[source_id]

    def source(self, source_id):
        return self.get_source(source_id)

    assessments = candidate_results
    get_candidate_results = candidate_results
    candidate_statuses = candidate_results
    current_candidate_results = candidate_results
    all_source_results = candidate_results
    list_candidates = candidate_results

    # ---- 待执行指令 ------------------------------------------------------

    def pending_instructions(self):
        items = [
            _public_instruction(instruction)
            for instruction in self._state["pending"]
        ]
        return IndexedList(items, "token", "source_id", "source")

    def pending_map(self, key="token"):
        return {item[key]: item for item in self.pending_instructions()}

    def instructions(self):
        return self.pending_instructions()

    def pending(self):
        return self.pending_instructions()

    def active_instructions(self):
        return self.pending_instructions()

    def list_instructions(self):
        return self.pending_instructions()

    def get_instruction(self, token):
        return self.pending_instructions()[token]

    pending_commands = pending_instructions
    get_pending_instructions = pending_instructions
    commands = pending_instructions
    command_list = pending_instructions
    current_instructions = pending_instructions
    list_pending = pending_instructions

    # ---- 已执行记录 ------------------------------------------------------

    def executed_records(self):
        items = [
            _public_execution(execution)
            for execution in self._state["executions"]
        ]
        return IndexedList(
            items,
            "token",
            "record_id",
            "id",
            "shot_id",
            "capture_id",
            "source_id",
            "source",
        )

    def execution_map(self, key="token"):
        return {item[key]: item for item in self.executed_records()}

    def executions(self):
        return self.executed_records()

    def records(self):
        return self.executed_records()

    def list_executions(self):
        return self.executed_records()

    execution_records = executed_records
    get_executed_records = executed_records
    executed_instructions = executed_records
    all_executions = executed_records
    executed = executed_records

    # ---- 组合读取 --------------------------------------------------------

    def read(self):
        candidates = self.candidate_results()
        pending = self.pending_instructions()
        executed = self.executed_records()
        current = [row for row in candidates if row["candidate"]]
        capacity = self.capacity_status()
        return {
            "candidates": candidates,
            "source_results": candidates,
            "current_candidates": current,
            "instructions": pending,
            "pending_instructions": pending,
            "pending": pending,
            "executions": executed,
            "executed_records": executed,
            **capacity,
        }

    def snapshot(self):
        return self.read()

    get_state = read
    get_snapshot = snapshot
    workbook_view = read

    # ---- 整批质量更正 ----------------------------------------------------

    def update_quality_flags(self, corrections=None, **kwargs):
        changes = _normalize_corrections(corrections, kwargs)

        def operation(state):
            measurement_ids = {item["id"] for item in state["measurements"]}
            # 先整体验证，再修改副本。
            for mid, _usable in changes:
                if mid not in measurement_ids:
                    raise FollowUpError(f"质量更正包含未知测量 ID：{mid}")

            by_id = {item["id"]: item for item in state["measurements"]}
            for mid, usable in changes:
                by_id[mid]["usable"] = usable

            assessment = _assess_state(state)
            plan = _reconcile_plan(state, assessment)
            return plan

        raw = self._commit(operation)
        revoked = IndexedList(
            [_public_instruction(item) for item in raw["revoked"]],
            "token",
            "source_id",
            "source",
        )
        issued = IndexedList(
            [_public_instruction(item) for item in raw["issued"]],
            "token",
            "source_id",
            "source",
        )
        pending = self.pending_instructions()
        return CorrectionResult(
            {
                "revoked": revoked,
                "issued": issued,
                # 明确的命名别名：
                "new": issued,
                "new_instructions": issued,
                "created": issued,
                "added": issued,
                "issued_instructions": issued,
                "revoked_instructions": revoked,
                "cancelled": revoked,
                "canceled": revoked,
                "invalidated": revoked,
                "instructions": pending,
                "pending_instructions": pending,
            }
        )

    # 常见方法名别名。
    apply_quality_updates = update_quality_flags
    apply_quality_update = update_quality_flags
    update_quality = update_quality_flags
    batch_update_quality = update_quality_flags
    apply_corrections = update_quality_flags
    apply_quality_corrections = update_quality_flags
    update_quality_markers = update_quality_flags
    apply_quality_flags = update_quality_flags
    set_quality_flags = update_quality_flags
    bulk_update_usable = update_quality_flags
    batch_correct_quality = update_quality_flags
    correct = update_quality_flags
    correct_quality_flags = update_quality_flags
    set_usable_flags = update_quality_flags
    update_flags = update_quality_flags
    update_usable = update_quality_flags
    update_usable_flags = update_quality_flags

    # ---- 登记执行 --------------------------------------------------------

    def record_execution(
        self,
        token=None,
        record_id=None,
        *extra_args,
        **kwargs,
    ):
        # 允许用 instruction= / instruction_token= 等关键字传令牌。
        if token is None:
            for key in ("instruction", "instruction_token", "command"):
                if key in kwargs:
                    token = kwargs.pop(key)
                    break

        object_record_id = None
        # 支持 record_execution({...}) / record_execution(instruction, shot_id)。
        if isinstance(token, dict):
            obj = token
            token = obj.get("token")
            object_record_id = (
                obj.get("record_id")
                or obj.get("shot_id")
                or obj.get("capture_id")
                or obj.get("实拍记录")
            )

        if extra_args:
            if record_id is None:
                record_id = extra_args[0]
            if len(extra_args) > 1:
                raise FollowUpError("登记执行只接受令牌和实拍记录 ID")

        id_values = [record_id, object_record_id]
        for key in (
            "shot_id",
            "capture_id",
            "observation_id",
            "exposure_id",
            "shot_record_id",
            "real_shot_id",
            "capture_record_id",
        ):
            if key in kwargs:
                id_values.append(kwargs.pop(key))

        normalized_ids = []
        for value in id_values:
            if value is not None:
                normalized_ids.append(value)
        if any(value != normalized_ids[0] for value in normalized_ids[1:]):
            raise FollowUpConflict("同一登记的实拍记录 ID 别名给出了冲突值")
        if normalized_ids:
            record_id = normalized_ids[0]

        if kwargs:
            raise FollowUpError(f"登记执行收到未定义参数：{sorted(kwargs)}")

        token = ident(token)
        record_id = ident(record_id)

        def operation(state):
            for execution in state["executions"]:
                if execution["token"] == token:
                    if execution["record_id"] == record_id:
                        # 幂等重试：原样返回已保存记录，绝不再次扣容量。
                        return {"execution": execution, "idempotent": True}
                    raise FollowUpConflict(
                        f"令牌 {token} 已对应另一实拍记录，不能改写"
                    )

            pending_by_token = {
                instruction["token"]: instruction
                for instruction in state["pending"]
            }
            if token not in pending_by_token:
                raise FollowUpConflict(
                    f"令牌 {token} 不是当前待执行指令（可能已执行、已撤销或未知）"
                )

            for execution in state["executions"]:
                if execution["record_id"] == record_id:
                    raise FollowUpConflict(
                        f"实拍记录 ID {record_id} 已用于另一令牌"
                    )

            instruction = pending_by_token[token]
            state["pending"] = [
                item for item in state["pending"] if item["token"] != token
            ]
            execution = {
                "token": token,
                "record_id": record_id,
                "instruction": copy.deepcopy(instruction),
            }
            state["executions"].append(execution)
            return {"execution": execution, "idempotent": False}

        result = self._commit(operation)
        return _public_execution(result["execution"])

    register_execution = record_execution
    register_shot = record_execution
    record_shot = record_execution
    register_capture = record_execution
    record_capture = record_execution
    register_observation = record_execution
    record_followup = record_execution
    mark_executed = record_execution
    complete_instruction = record_execution
    execute_instruction = record_execution


# ---------------------------------------------------------------------------
# 模块级便捷函数和名称别名
# ---------------------------------------------------------------------------

Workbook = FollowUpWorkbook
WorkBook = FollowUpWorkbook
FollowUpWorkBook = FollowUpWorkbook
FollowupWorkBook = FollowUpWorkbook
FollowupWorkbook = FollowUpWorkbook
OfflineFollowUpWorkbook = FollowUpWorkbook
FollowUpStore = FollowUpWorkbook
FollowupStore = FollowUpWorkbook


def create_workbook(path, data=None, measurements=None, **options):
    return FollowUpWorkbook.create(path, data, measurements, **options)


def initialize_workbook(path, data=None, measurements=None, **options):
    return FollowUpWorkbook.initialize(path, data, measurements, **options)


def init_workbook(path, data=None, measurements=None, **options):
    return FollowUpWorkbook.init(path, data, measurements, **options)


def open_workbook(path):
    return FollowUpWorkbook.open(path)


def load_workbook(path):
    return FollowUpWorkbook.load(path)


def followup_workbook(path=None, data=None, measurements=None, **options):
    """无状态风格入口：有批次则创建，无批次则重开。"""
    if data is not None or measurements is not None or options:
        return create_workbook(path, data, measurements, **options)
    return open_workbook(path)


create_followup_workbook = create_workbook
create_followup = create_workbook
initialize_followup_workbook = initialize_workbook
open_followup_workbook = open_workbook
open_followup = open_workbook
load_followup_workbook = load_workbook
load_followup = load_workbook
FollowUp = FollowUpWorkbook
FollowUpBook = FollowUpWorkbook
PhotometryWorkbook = FollowUpWorkbook
