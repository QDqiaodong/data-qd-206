"""测光发布件与变暗判据的离线验收测试(仅标准库)。

覆盖:
* 任务书指定的 A/B/C/D 四帧主剧本、+0.5 变体、A 实测/C 极限变体;
* 无证明/证明被几何候选阻止(可区分、可追溯检测 ID);
* 星等可破平局但几何并列时不签发实测;
* 无基准连通分量可比较、无共同标尺不可比较、双定值分量可比较;
* 需排除两条转移的矛盾分量不阻塞别组;无矛盾分量只有“不排除”假设;
* 坏引用、重复身份、逆时间、重排不变、输入不被修改;
* 发布件经 JSON 往返后在独立进程里得到同样结论;
* 发布件不携带原始曝光行、坐标、候选边表或整份 photometry。
"""
import copy
import json
import subprocess
import sys
import unittest
from pathlib import Path

from star_crossmatch import (
    PHOTOMETRIC_RELEASE_FORMAT,
    PhotometricRelease,
    assess_fade,
    build_photometric_release,
    process_exposures,
)

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def ref(rid, ra, dec, mag, cal=False):
    return {"id": rid, "ra": ra, "dec": dec, "magnitude": mag, "calibrator": cal}


def obs(oid, ra, dec, sigma, instr):
    return {
        "id": oid,
        "ra": ra,
        "dec": dec,
        "sigma": sigma,
        "instrumental_magnitude": instr,
    }


def transfer(tid, src, dst, delta):
    return {"id": tid, "from_exposure": src, "to_exposure": dst, "delta": delta}


def limit_row(eid, rid, bound):
    return {
        "exposure": eid,
        "reference": rid,
        "limiting_instrumental_magnitude": bound,
    }


def make_batch(transfers, limits, *, cd_delta="1", d_extra=None, c_obs=None):
    """任务书 A/B/C/D 主剧本的可复用构造器。

    A 含三个等偏移(+2)定标锚点;S 是目标变源。C 帧默认含 S 的几何确定
    实测(仪器星等 16);D 帧默认不含 S 的几何候选。
    """
    references = [
        ref("C1", "10", "0", "15", True),
        ref("C2", "20", "0", "15", True),
        ref("C3", "30", "0", "15", True),
        ref("S", "50", "0", "18"),
        ref("T", "50.0002", "0", "18"),
    ]
    d_observations = [obs("d1", "70", "0", "0.001", "13")]
    if d_extra:
        d_observations.extend(d_extra)
    c_observations = c_obs if c_obs is not None else [
        obs("sC", "50", "0", "0.001", "16")
    ]
    photometry = {
        "reference": references,
        "threshold": "1",
        "anchor_tolerance": "0.5",
        "exposures": [
            {"id": "A", "observations": [
                obs("k1", "10", "0", "0.001", "13"),
                obs("k2", "20", "0", "0.001", "13"),
                obs("k3", "30", "0", "0.001", "13"),
            ]},
            {"id": "B", "observations": [obs("b1", "80", "0", "0.001", "13")]},
            {"id": "C", "observations": c_observations},
            {"id": "D", "observations": d_observations},
        ],
    }
    epochs = [
        {"exposure": "A", "mjd": "1"},
        {"exposure": "B", "mjd": "2"},
        {"exposure": "C", "mjd": "3"},
        {"exposure": "D", "mjd": "4"},
    ]
    if transfers is None:
        transfers = [
            transfer("t1", "A", "B", "1"),
            transfer("t2", "B", "C", "1"),
            transfer("t3", "A", "C", "3"),
            transfer("t4", "C", "D", cd_delta),
        ]
    return {
        "photometry": photometry,
        "epochs": epochs,
        "transfers": transfers,
        "limits": limits,
    }


def fade_query(bright, later, minimum, reference="S"):
    return {
        "reference": reference,
        "bright_exposure": bright,
        "later_exposure": later,
        "minimum_fade": minimum,
    }


class MandatoryScenario(unittest.TestCase):
    def setUp(self):
        self.data = make_batch(
            None, [limit_row("D", "S", "16")]
        )
        self.release = build_photometric_release(self.data)

    def test_three_single_exclusion_hypotheses_and_unpublished_zeros(self):
        comp = self.release["components"][0]
        self.assertEqual(comp["status"], "solvable")
        self.assertEqual(
            [h["excluded_transfer"] for h in comp["hypotheses"]],
            ["t1", "t2", "t3"],
        )
        statuses = {row["exposure"]: row for row in self.release["exposures"]}
        self.assertEqual(statuses["A"]["zero_point_status"], "published")
        self.assertEqual(statuses["A"]["zero_point"], "2")
        self.assertTrue(statuses["A"]["is_anchor"])
        for eid in ("B", "C", "D"):
            self.assertEqual(statuses[eid]["zero_point_status"], "indeterminate")
            self.assertIsNone(statuses[eid]["zero_point"])
            self.assertFalse(statuses[eid]["is_anchor"])

    def test_measurement_and_clean_proof(self):
        key = ("C", "S")
        measured = {
            (m["exposure"], m["reference"]): m
            for m in self.release["measurements"]
        }
        self.assertEqual(
            measured[key]["observation"], "sC"
        )
        self.assertEqual(measured[key]["instrumental_magnitude"], "16")
        proofs = {
            (p["exposure"], p["reference"]): p
            for p in self.release["sensitivity_proofs"]
        }
        self.assertTrue(proofs[("D", "S")]["usable"])
        self.assertEqual(proofs[("D", "S")]["blocking_detections"], [])

    def test_fade_confirmed_with_bound_exactly_one(self):
        result = assess_fade(self.release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "fade_confirmed")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "1")
        self.assertEqual(
            result["hypothesis_evidence"],
            [{"component": "component-1", "excluded_transfer": "t1"}],
        )

    def test_half_delta_variant_lowers_bound_to_half(self):
        data = make_batch(None, [limit_row("D", "S", "16")], cd_delta="0.5")
        release = build_photometric_release(data)
        self.assertEqual(
            [h["excluded_transfer"] for h in release["components"][0]["hypotheses"]],
            ["t1", "t2", "t3"],
        )
        statuses = {row["exposure"]: row for row in release["exposures"]}
        for eid in ("B", "C", "D"):
            self.assertEqual(statuses[eid]["zero_point_status"], "indeterminate")
        result = assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "insufficient_evidence")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "0.5")

    def test_a_measurement_c_limit_confirms_at_two(self):
        # 独立构造:A 帧含 S 的实测(16),C 帧没有 S 的几何候选但有清洁证明。
        data = make_batch(
            None,
            [limit_row("C", "S", "16")],
            c_obs=[obs("c1", "80", "0", "0.001", "13")],
        )
        data["photometry"]["exposures"][0]["observations"].append(
            obs("sA", "50", "0", "0.001", "16")
        )
        release = build_photometric_release(data)
        result = assess_fade(release, fade_query("A", "C", "2"))
        self.assertEqual(result["status"], "fade_confirmed")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "2")
        # 三个假设给出的差值分别为 3、3、2,最小者来自排除 t3。
        self.assertEqual(
            result["hypothesis_evidence"],
            [{"component": "component-1", "excluded_transfer": "t3"}],
        )
        # 更严的幅度 2.0001 不应确认,但下界仍是 2。
        weak = assess_fade(release, fade_query("A", "C", "2.0001"))
        self.assertEqual(weak["status"], "insufficient_evidence")
        self.assertEqual(weak["guaranteed_fade_lower_bound"], "2")


class EvidenceRejections(unittest.TestCase):
    def test_missing_sensitivity_proof(self):
        data = make_batch(None, [])
        release = build_photometric_release(data)
        result = assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "no_sensitivity_proof")
        self.assertEqual(result["detections"], [])

    def test_blocked_proof_names_detection_and_is_distinct_status(self):
        # dS 位于 T 的位置:它同时是 S 的原始几何候选(角距 0.0002 度),
        # 但几何最优方案始终把它配给 T(代价 0 优于配给 S 的 720 毫角秒),
        # 从未成为 S 的最终配对,却仍永久阻止 S 的未检出证明。
        blocking_dets = [obs("dS", "50.0002", "0", "0.001", "13")]
        data = make_batch(
            None,
            [limit_row("D", "S", "16")],
            d_extra=blocking_dets,
        )
        release = build_photometric_release(data)
        proofs = {
            (p["exposure"], p["reference"]): p
            for p in release["sensitivity_proofs"]
        }
        self.assertFalse(proofs[("D", "S")]["usable"])
        self.assertEqual(proofs[("D", "S")]["blocking_detections"], ["dS"])
        result = assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "sensitivity_proof_blocked")
        self.assertEqual(result["detections"], ["dS"])

    def test_no_trusted_measurement_traces_candidate_detections(self):
        # C 帧两条检测在 S 处几何并列,没有任何一条在全部几何最优方案中锁定 S。
        tied = [
            obs("p1", "50", "0", "0.001", "16"),
            obs("p2", "50", "0", "0.001", "17"),
        ]
        data = make_batch(
            None, [limit_row("D", "S", "16")], c_obs=tied
        )
        release = build_photometric_release(data)
        self.assertNotIn(
            ("C", "S"),
            {(m["exposure"], m["reference"]) for m in release["measurements"]},
        )
        result = assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "no_trusted_measurement")
        self.assertEqual(result["detections"], ["p1", "p2"])
        # 同样信息必须在 JSON 往返后的独立发布件里保留。
        rebuilt = json.loads(json.dumps(release))
        result2 = assess_fade(rebuilt, fade_query("C", "D", "1"))
        self.assertEqual(result2, result)

    def test_magnitude_tie_break_in_old_api_does_not_create_measurement(self):
        # 独立可定标帧:三颗定标星给有效零点;N1/N2 同位、p1/p2 同几何代价。
        references = [
            ref("K1", "10", "0", "15", True),
            ref("K2", "20", "0", "15", True),
            ref("K3", "30", "0", "15", True),
            ref("N1", "60", "0", "15"),
            ref("N2", "60", "0", "17"),
        ]
        photometry = {
            "reference": references,
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "E", "observations": [
                    obs("k1", "10", "0", "0.001", "13"),
                    obs("k2", "20", "0", "0.001", "13"),
                    obs("k3", "30", "0", "0.001", "13"),
                    obs("p1", "60", "0", "0.001", "15"),
                    obs("p2", "60", "0", "0.001", "13"),
                ]}
            ],
        }
        legacy = process_exposures(photometry)["exposures"][0]
        detmap = {d["observation"]: d for d in legacy["detections"]}
        # 原 API 在有效零点下按星等残差破平局:p2(校正 15)配 N1、p1(校正 17)配 N2。
        self.assertEqual(legacy["calibration"]["zero_point"], "2")
        self.assertEqual(detmap["p1"]["selected"], "N2")
        self.assertEqual(detmap["p2"]["selected"], "N1")
        data = {
            "photometry": photometry,
            "epochs": [{"exposure": "E", "mjd": "1"}],
            "transfers": [],
            "limits": [],
        }
        release = build_photometric_release(data)
        # 三颗定标星的强制配对仍是实测;N1/N2 几何并列,不得签发实测。
        measured_refs = {m["reference"] for m in release["measurements"]}
        self.assertEqual(measured_refs, {"K1", "K2", "K3"})
        identities = {
            (g["exposure"], g["reference"]): g["candidate_detections"]
            for g in release["geometric_identities"]
        }
        self.assertEqual(identities[("E", "N1")], ["p1", "p2"])
        self.assertEqual(identities[("E", "N2")], ["p1", "p2"])


class CalibrationComponents(unittest.TestCase):
    def test_consistent_component_keeps_only_no_exclusion_hypothesis(self):
        # 自洽三角形:三条转移互相印证,只有 null 排除一种假设。
        transfers = [
            transfer("u1", "A", "B", "1"),
            transfer("u2", "B", "C", "1"),
            transfer("u3", "A", "C", "2"),
        ]
        data = make_batch(transfers, [])
        release = build_photometric_release(data)
        comp = release["components"][0]
        self.assertEqual(
            [h["excluded_transfer"] for h in comp["hypotheses"]], [None]
        )
        # A 是锚点,全链绝对值都可发布。
        values = {row["exposure"]: row["zero_point"] for row in release["exposures"]}
        self.assertEqual(values, {"A": "2", "B": "3", "C": "4", "D": None})
        self.assertEqual(
            next(row for row in release["exposures"] if row["exposure"] == "D")[
                "zero_point_status"
            ],
            "unanchored",
        )

    def test_unanchored_connected_frames_remain_comparable(self):
        # E、F 无锚点但由一条自洽转移连通;差值确定,可直接作差。
        photometry = {
            "reference": [ref("S", "10", "0", "18")],
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "E", "observations": [obs("sE", "10", "0", "0.001", "16")]},
                {"id": "F", "observations": [obs("f1", "80", "0", "0.001", "13")]},
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "1"},
                {"exposure": "F", "mjd": "2"},
            ],
            "transfers": [transfer("x1", "E", "F", "0.25")],
            "limits": [limit_row("F", "S", "16")],
        }
        release = build_photometric_release(data)
        for row in release["exposures"]:
            self.assertEqual(row["zero_point_status"], "unanchored")
            self.assertIsNone(row["zero_point"])
        result = assess_fade(release, fade_query("E", "F", "1"))
        # 16 - 16 + 0.25 = 0.25。
        self.assertEqual(result["status"], "insufficient_evidence")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "0.25")

    def test_two_anchored_components_are_comparable_without_link(self):
        photometry = {
            "reference": [
                ref("K1", "10", "0", "15", True),
                ref("K2", "20", "0", "15", True),
                ref("K3", "30", "0", "15", True),
                ref("S", "50", "0", "18"),
            ],
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "E", "observations": [
                    obs("k1", "10", "0", "0.001", "13"),
                    obs("k2", "20", "0", "0.001", "13"),
                    obs("k3", "30", "0", "0.001", "13"),
                    obs("sE", "50", "0", "0.001", "16"),
                ]},
                {"id": "F", "observations": [
                    obs("j1", "10", "0", "0.001", "12"),
                    obs("j2", "20", "0", "0.001", "12"),
                    obs("j3", "30", "0", "0.001", "12"),
                    obs("f1", "80", "0", "0.001", "13"),
                ]},
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "1"},
                {"exposure": "F", "mjd": "2"},
            ],
            "transfers": [],
            "limits": [limit_row("F", "S", "16")],
        }
        release = build_photometric_release(data)
        # 两个独立分量,各自有确定绝对零点 +2 与 +3。
        self.assertEqual(len(release["components"]), 2)
        values = {row["exposure"]: row["zero_point"] for row in release["exposures"]}
        self.assertEqual(values, {"E": "2", "F": "3"})
        result = assess_fade(release, fade_query("E", "F", "2"))
        # 16 - 16 + (3 - 2) = 1,达不到 2。
        self.assertEqual(result["status"], "insufficient_evidence")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "1")
        confirmed = assess_fade(release, fade_query("E", "F", "1"))
        self.assertEqual(confirmed["status"], "fade_confirmed")
        self.assertEqual(
            confirmed["hypothesis_evidence"],
            [
                {"component": "component-1", "excluded_transfer": None},
                {"component": "component-2", "excluded_transfer": None},
            ],
        )

    def test_no_common_scale_is_not_comparable(self):
        # E、F 分属两个无锚分量,无转移相连:零点差不受约束。
        photometry = {
            "reference": [ref("S", "10", "0", "18")],
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "E", "observations": [obs("sE", "10", "0", "0.001", "16")]},
                {"id": "F", "observations": [obs("f1", "80", "0", "0.001", "13")]},
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "1"},
                {"exposure": "F", "mjd": "2"},
            ],
            "transfers": [],
            "limits": [limit_row("F", "S", "16")],
        }
        release = build_photometric_release(data)
        result = assess_fade(release, fade_query("E", "F", "1"))
        self.assertEqual(result["status"], "not_comparable")
        self.assertNotIn("guaranteed_fade_lower_bound", result)
        self.assertEqual(
            result["hypothesis_evidence"],
            [
                {"component": "component-1", "excluded_transfer": None},
                {"component": "component-2", "excluded_transfer": None},
            ],
        )

    def test_two_removal_component_does_not_block_other_component(self):
        # 分量 X:锚点 A,四条边含两个独立矛盾,删任何一条都仍有矛盾 => 无解。
        good = [
            transfer("g1", "E", "F", "1"),
            transfer("g2", "F", "G", "1"),
            transfer("g3", "E", "G", "2"),
        ]
        bad = [
            transfer("b1", "A", "B", "1"),
            transfer("b2", "B", "C", "1"),
            transfer("b3", "A", "C", "3"),  # 与 b1+b2=2 矛盾
            transfer("b4", "A", "D", "1"),
            transfer("b5", "D", "B", "5"),  # 与 A->B=1 矛盾
        ]
        references = [
            ref("K1", "10", "0", "15", True),
            ref("K2", "20", "0", "15", True),
            ref("K3", "30", "0", "15", True),
            ref("S", "50", "0", "18"),
        ]
        far = "80"
        photometry = {
            "reference": references,
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "A", "observations": [
                    obs("k1", "10", "0", "0.001", "13"),
                    obs("k2", "20", "0", "0.001", "13"),
                    obs("k3", "30", "0", "0.001", "13"),
                ]},
                {"id": "B", "observations": [
                    obs("b1o", far, "0", "0.001", "13"),
                    obs("sB", "50", "0", "0.001", "16"),
                ]},
                {"id": "C", "observations": [obs("c1o", far, "0", "0.001", "13")]},
                {"id": "D", "observations": [obs("d1o", far, "0", "0.001", "13")]},
                {"id": "E", "observations": [obs("sE", "50", "0", "0.001", "16")]},
                {"id": "F", "observations": [obs("f1o", far, "0", "0.001", "13")]},
                {"id": "G", "observations": [obs("g1o", far, "0", "0.001", "13")]},
            ],
        }
        epochs = [
            {"exposure": eid, "mjd": str(index + 1)}
            for index, eid in enumerate("ABCDEFG")
        ]
        data = {
            "photometry": photometry,
            "epochs": epochs,
            "transfers": good + bad,
            "limits": [limit_row("F", "S", "16")],
        }
        release = build_photometric_release(data)
        comps = {c["component"]: c for c in release["components"]}
        bad_comp = next(c for c in release["components"] if "A" in c["exposures"])
        good_comp = next(c for c in release["components"] if "E" in c["exposures"])
        self.assertEqual(bad_comp["status"], "unsolvable")
        self.assertEqual(bad_comp["hypotheses"], [])
        self.assertEqual(good_comp["status"], "solvable")
        # 无解散分量内帧(含锚点 A)状态为 unsolvable 且不发布零点。
        for eid in "ABCD":
            row = next(r for r in release["exposures"] if r["exposure"] == eid)
            self.assertEqual(row["zero_point_status"], "unsolvable")
            self.assertIsNone(row["zero_point"])
        # 好分量照常发布并支持查询:16-16+1 = 1,确认下界 1。
        result = assess_fade(release, fade_query("E", "F", "1"))
        self.assertEqual(result["status"], "fade_confirmed")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "1")
        # 触碰到无解散分量的查询返回定标无解,且不抹掉别组成果。
        dead = assess_fade(
            release,
            {
                "reference": "S",
                "bright_exposure": "B",
                "later_exposure": "F",
                "minimum_fade": "1",
            },
        )
        self.assertEqual(dead["status"], "calibration_unsolvable")

    def test_parallel_transfers_are_independent_debug_units(self):
        # 两条独立的 A->B 记录互相矛盾(不得先合并):分别排除各自可行,
        # 其余约束自洽;B 与 A、C 之间没有别的路径。
        transfers = [
            transfer("p1", "A", "B", "1"),
            transfer("p2", "A", "B", "2"),
            transfer("p3", "A", "C", "3"),
        ]
        data = make_batch(transfers, [])
        release = build_photometric_release(data)
        self.assertEqual(
            sorted(h["excluded_transfer"] for h in release["components"][0]["hypotheses"]),
            ["p1", "p2"],
        )

    def test_disconnected_split_with_two_anchors_still_comparable(self):
        # 同一原分量含两个绝对基准 E、F,桥接转移与二者都冲突:唯一可行假设是
        # 删桥,假设内两帧不互通,但两侧都有确定绝对值,仍可按绝对值作差。
        anchors = [
            ref("K1", "10", "0", "15", True),
            ref("K2", "20", "0", "15", True),
            ref("K3", "30", "0", "15", True),
            ref("L1", "11", "0", "15", True),
            ref("L2", "21", "0", "15", True),
            ref("L3", "31", "0", "15", True),
            ref("S", "50", "0", "18"),
        ]
        photometry = {
            "reference": anchors,
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {"id": "E", "observations": [
                    obs("k1", "10", "0", "0.001", "13"),
                    obs("k2", "20", "0", "0.001", "13"),
                    obs("k3", "30", "0", "0.001", "13"),
                    obs("sE", "50", "0", "0.001", "16"),
                ]},
                {"id": "F", "observations": [
                    obs("l1", "11", "0", "0.001", "12"),
                    obs("l2", "21", "0", "0.001", "12"),
                    obs("l3", "31", "0", "0.001", "12"),
                    obs("f1", "80", "0", "0.001", "13"),
                ]},
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "1"},
                {"exposure": "F", "mjd": "2"},
            ],
            # 真实零点差为 3-2=1,桥却声称 9,与两个锚点都冲突。
            "transfers": [transfer("q1", "E", "F", "9")],
            "limits": [limit_row("F", "S", "16")],
        }
        release = build_photometric_release(data)
        comp = release["components"][0]
        self.assertEqual(comp["status"], "solvable")
        self.assertEqual(
            [h["excluded_transfer"] for h in comp["hypotheses"]], ["q1"]
        )
        values = {row["exposure"]: row["zero_point"] for row in release["exposures"]}
        self.assertEqual(values, {"E": "2", "F": "3"})
        # 16 - 16 + (3 - 2) = 1。
        result = assess_fade(release, fade_query("E", "F", "1"))
        self.assertEqual(result["status"], "fade_confirmed")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "1")
        self.assertEqual(
            result["hypothesis_evidence"],
            [{"component": "component-1", "excluded_transfer": "q1"}],
        )


class FadeCalibrationGate(unittest.TestCase):
    def test_unsolvable_component_short_circuits_before_bound(self):
        data = make_batch(
            [
                transfer("b1", "A", "B", "1"),
                transfer("b2", "B", "C", "1"),
                transfer("b3", "A", "C", "3"),
                transfer("b4", "A", "D", "1"),
                transfer("b5", "D", "B", "5"),
            ],
            [limit_row("D", "S", "16")],
        )
        release = build_photometric_release(data)
        result = assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "calibration_unsolvable")
        self.assertNotIn("guaranteed_fade_lower_bound", result)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.data = make_batch(None, [limit_row("D", "S", "16")])

    def expect_value_error(self, mutate):
        data = copy.deepcopy(self.data)
        mutate(data)
        with self.assertRaises(ValueError):
            build_photometric_release(data)

    def test_top_level_fields_strict(self):
        self.expect_value_error(lambda d: d.pop("limits"))
        self.expect_value_error(lambda d: d.update(extra=1))

    def test_new_size_caps(self):
        def too_many_exposures(d):
            base = d["photometry"]["exposures"][0]
            d["photometry"]["exposures"] = [
                {"id": f"X{i}", "observations": copy.deepcopy(base["observations"])}
                for i in range(33)
            ]
            d["epochs"] = [
                {"exposure": f"X{i}", "mjd": "1"} for i in range(33)
            ]
        self.expect_value_error(too_many_exposures)

        def too_many_refs(d):
            refs = d["photometry"]["reference"]
            template = refs[0].copy()
            for i in range(257 - len(refs)):
                row = template.copy()
                row["id"] = f"R{i}"
                refs.append(row)
        self.expect_value_error(too_many_refs)

        def too_many_detections(d):
            d["photometry"]["reference"] = [
                ref("S", "10", "0", "18")
            ]
            d["photometry"]["exposures"] = [
                {
                    "id": "A",
                    "observations": [
                        obs(f"o{i}", "10", "0", "0.1", "13") for i in range(2001)
                    ],
                }
            ]
            d["epochs"] = [{"exposure": "A", "mjd": "1"}]
        self.expect_value_error(too_many_detections)

    def test_transfer_rules(self):
        self.expect_value_error(lambda d: d["transfers"].append(
            transfer("t1", "A", "B", "9")))  # 重复 ID
        self.expect_value_error(lambda d: d["transfers"].append(
            transfer("t9", "A", "A", "1")))  # 同帧
        self.expect_value_error(lambda d: d["transfers"].append(
            transfer("t9", "A", "Z", "1")))  # 未知帧
        self.expect_value_error(lambda d: d["transfers"][0].update(delta=1))  # 类型

    def test_epoch_rules(self):
        self.expect_value_error(lambda d: d["epochs"].append(
            {"exposure": "A", "mjd": "5"}))  # 每帧恰一行
        self.expect_value_error(lambda d: d["epochs"].append(
            {"exposure": "Z", "mjd": "5"}))  # 未知帧
        self.expect_value_error(lambda d: d["epochs"][0].update(mjd="-1"))
        self.expect_value_error(lambda d: d["epochs"][0].pop("mjd"))

    def test_limit_rules(self):
        self.expect_value_error(lambda d: d["limits"].append(
            limit_row("D", "S", "17")))  # 组合重复
        self.expect_value_error(lambda d: d["limits"].append(
            limit_row("Z", "S", "17")))
        self.expect_value_error(lambda d: d["limits"].append(
            limit_row("D", "Z", "17")))
        self.expect_value_error(
            lambda d: d["limits"][0].update(limiting_instrumental_magnitude=17.0))

    def test_query_validation(self):
        release = build_photometric_release(self.data)
        with self.assertRaises(ValueError):
            assess_fade(release, fade_query("C", "D", "0"))  # 幅度必须为正
        with self.assertRaises(ValueError):
            assess_fade(release, fade_query("D", "C", "1"))  # 逆时间
        with self.assertRaises(ValueError):  # 同刻也不允许(必须严格小于)
            assess_fade(
                release,
                {
                    "reference": "S",
                    "bright_exposure": "A",
                    "later_exposure": "A",
                    "minimum_fade": "1",
                },
            )
        with self.assertRaises(ValueError):
            assess_fade(release, fade_query("C", "Z", "1"))
        with self.assertRaises(ValueError):
            assess_fade(release, fade_query("Z", "D", "1"))
        with self.assertRaises(ValueError):
            assess_fade(release, {**fade_query("C", "D", "1"), "extra": 1})

    def test_unsupported_release_format(self):
        with self.assertRaises(ValueError):
            assess_fade({"format": "something-else/9"}, fade_query("C", "D", "1"))
        with self.assertRaises(ValueError):
            assess_fade("not-a-dict", fade_query("C", "D", "1"))
        payload = json.loads(json.dumps(build_photometric_release(self.data)))
        payload["format"] = PHOTOMETRIC_RELEASE_FORMAT + "-tampered"
        with self.assertRaises(ValueError):
            assess_fade(payload, fade_query("C", "D", "1"))

    def test_input_not_mutated(self):
        snapshot = copy.deepcopy(self.data)
        build_photometric_release(self.data)
        self.assertEqual(self.data, snapshot)
        release = build_photometric_release(self.data)
        payload_snapshot = copy.deepcopy(dict(release))
        assess_fade(release, fade_query("C", "D", "1"))
        self.assertEqual(dict(release), payload_snapshot)


class OrderInvariance(unittest.TestCase):
    def test_reordering_inputs_does_not_change_release_or_result(self):
        data = make_batch(None, [limit_row("D", "S", "16")])
        shuffled = copy.deepcopy(data)
        shuffled["photometry"]["reference"] = list(
            reversed(shuffled["photometry"]["reference"])
        )
        shuffled["photometry"]["exposures"] = list(
            reversed(shuffled["photometry"]["exposures"])
        )
        for exp in shuffled["photometry"]["exposures"]:
            exp["observations"] = list(reversed(exp["observations"]))
        shuffled["transfers"] = list(reversed(shuffled["transfers"]))
        shuffled["limits"] = list(reversed(shuffled["limits"]))
        shuffled["epochs"] = list(reversed(shuffled["epochs"]))
        r1 = build_photometric_release(data)
        r2 = build_photometric_release(shuffled)
        self.assertEqual(json.loads(json.dumps(r1)), json.loads(json.dumps(r2)))
        query = fade_query("C", "D", "1")
        self.assertEqual(assess_fade(r1, query), assess_fade(r2, query))


class OfflineHandover(unittest.TestCase):
    CHILD = r"""
import json, sys
sys.path.insert(0, sys.argv[1])
from star_crossmatch import assess_fade
release = json.loads(sys.stdin.read())
query = json.loads(sys.argv[2])
print(json.dumps(assess_fade(release, query), ensure_ascii=False, sort_keys=True))
"""

    def run_child(self, release, query):
        payload = json.dumps(release, ensure_ascii=False, sort_keys=True)
        proc = subprocess.run(
            [sys.executable, "-I", "-c", self.CHILD, str(SRC), json.dumps(query)],
            input=payload,
            capture_output=True,
            text=True,
            check=True,
        )
        return json.loads(proc.stdout)

    def test_json_roundtrip_in_independent_process(self):
        data = make_batch(None, [limit_row("D", "S", "16")])
        release = build_photometric_release(data)
        query = fade_query("C", "D", "1")
        in_process = assess_fade(json.loads(json.dumps(release)), query)
        in_new_process = self.run_child(
            json.loads(json.dumps(release)), query
        )
        self.assertEqual(in_new_process, in_process)
        self.assertEqual(in_new_process["status"], "fade_confirmed")
        self.assertEqual(in_new_process["guaranteed_fade_lower_bound"], "1")

    def test_half_variant_in_independent_process(self):
        data = make_batch(None, [limit_row("D", "S", "16")], cd_delta="0.5")
        release = json.loads(json.dumps(build_photometric_release(data)))
        result = self.run_child(release, fade_query("C", "D", "1"))
        self.assertEqual(result["status"], "insufficient_evidence")
        self.assertEqual(result["guaranteed_fade_lower_bound"], "0.5")


class ReleaseContract(unittest.TestCase):
    def test_versioned_json_payload_without_raw_observations(self):
        data = make_batch(None, [limit_row("D", "S", "16")])
        release = build_photometric_release(data)
        self.assertIsInstance(release, PhotometricRelease)
        self.assertEqual(release["format"], PHOTOMETRIC_RELEASE_FORMAT)
        text = release.to_json()
        again = json.loads(text)
        self.assertEqual(again["format"], PHOTOMETRIC_RELEASE_FORMAT)
        # 不得转交原始曝光行、坐标、匹配候选边表或整份 photometry。
        forbidden_keys = {
            "ra", "dec", "sigma", "observations", "photometry",
            "candidates", "distance", "cost", "threshold",
            "anchor_tolerance", "calibrator",
        }

        def collect_keys(node, found):
            if isinstance(node, dict):
                found.update(node.keys())
                for value in node.values():
                    collect_keys(value, found)
            elif isinstance(node, list):
                for value in node:
                    collect_keys(value, found)

        found_keys = set()
        collect_keys(again, found_keys)
        self.assertEqual(found_keys & forbidden_keys, set())
        # 所有数值字段均为十进制字符串(或 null),没有裸浮点。
        def walk(node):
            if isinstance(node, dict):
                for key, value in node.items():
                    self.assertNotIsInstance(value, float)
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    self.assertNotIsInstance(value, float)
                    walk(value)
        walk(again)
        # 假设身份稳定且带版本分量编号。
        identities = [
            h["hypothesis"] for h in again["components"][0]["hypotheses"]
        ]
        self.assertEqual(
            identities,
            ["component-1#h1", "component-1#h2", "component-1#h3"],
        )

    def test_tampered_payload_is_rejected_on_load(self):
        data = make_batch(None, [limit_row("D", "S", "16")])
        payload = json.loads(json.dumps(build_photometric_release(data)))
        # 人为制造不可行假设:发布件装载时复算约束必须发现。
        payload["components"][0]["hypotheses"][0]["excluded_transfer"] = None
        with self.assertRaises(ValueError):
            assess_fade(payload, fade_query("C", "D", "1"))
        # 假设引用不存在的转移也必须拒绝。
        payload2 = json.loads(json.dumps(build_photometric_release(data)))
        payload2["components"][0]["hypotheses"][0]["excluded_transfer"] = "nope"
        with self.assertRaises(ValueError):
            assess_fade(payload2, fade_query("C", "D", "1"))


if __name__ == "__main__":
    unittest.main()
