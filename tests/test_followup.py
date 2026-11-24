import copy
import json
import tempfile
import unittest
from pathlib import Path

from star_crossmatch.followup import (
    FollowUpConflict,
    create_workbook,
    open_workbook,
)


def acceptance_batch():
    priorities = {"X": 1, "B": 2, "D": 3, "A": 4, "E": 5}
    counts = {
        "X": (8, "bright"),
        "B": (5, "bright"),
        "D": (7, "bright"),
        "A": (6, "bright"),
        "E": (5, "faint"),
    }
    sources = [
        {
            "id": sid,
            "baseline_magnitude": "15",
            "priority": priority,
        }
        for sid, priority in priorities.items()
    ]
    measurements = []
    for sid, (count, kind) in counts.items():
        for index in range(count):
            measurements.append(
                {
                    "id": f"{sid.lower()}{index}",
                    "source_id": sid,
                    # 各源使用互不重复的时间区间。
                    "mjd": str(priorities[sid] * 100 + index),
                    "magnitude": "14.9" if kind == "bright" else "15.1",
                    "sigma": "0.01",
                    "usable": True,
                }
            )
    return {
        "sources": sources,
        "measurements": measurements,
        "min_valid_times": 5,
        "fdr_threshold": "0.04",
        "capacity": 2,
    }


class RequiredAcceptance(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "followup.json"
        self.wb = create_workbook(self.path, acceptance_batch())

    def test_full_close_reopen_correction_and_conflicts(self):
        candidates = self.wb.candidate_results()
        self.assertEqual(candidates["B"]["p"], "1/32")
        self.assertEqual(candidates["B"]["q"], "5/128")

        pending = self.wb.pending_instructions()
        self.assertEqual([item["source_id"] for item in pending], ["X", "B"])
        x_token = pending["X"]["token"]
        b_token = pending["B"]["token"]
        self.assertEqual(pending["X"]["q"], "5/256")
        self.assertEqual(pending["B"]["q"], "5/128")

        # 保留“原登记输入”和第一次回执。
        submission = (x_token, "shot-X-001")
        receipt = self.wb.record_execution(*submission)
        self.assertEqual(receipt["instruction"]["source_id"], "X")
        self.assertEqual(receipt["instruction"]["q"], "5/256")
        self.wb.close()

        wb = open_workbook(self.path)
        self.assertEqual(wb.execution_count, 1)

        # 已执行令牌不能换实拍 ID；同一实拍 ID 也不能用于 B 令牌。
        with self.assertRaises((ValueError, FollowUpConflict)):
            wb.record_execution(x_token, "shot-X-changed")
        with self.assertRaises((ValueError, FollowUpConflict)):
            wb.record_execution(b_token, "shot-X-001")

        # 屏蔽 A 的两个时刻：A 的非相等时刻 6 -> 4，变为证据不足。
        changed = wb.update_quality_flags(
            [
                {"measurement_id": "a0", "usable": False},
                {"measurement_id": "a1", "usable": False},
            ]
        )
        self.assertEqual(
            [item["source_id"] for item in changed["revoked"]], ["B"]
        )
        self.assertEqual(changed["revoked"][0]["token"], b_token)
        self.assertEqual(
            [item["source_id"] for item in changed["issued"]], ["D"]
        )
        d_token = changed["issued"][0]["token"]

        current = wb.candidate_results()
        self.assertEqual(current["A"]["n"], 4)
        self.assertEqual(current["A"]["k"], 4)
        self.assertIsNone(current["A"]["p"])
        self.assertIsNone(current["A"]["q"])
        self.assertEqual(current["B"]["n"], 5)
        self.assertEqual(current["B"]["k"], 5)
        self.assertEqual(current["B"]["p"], "1/32")
        self.assertEqual(current["B"]["q"], "1/24")
        self.assertEqual(current["X"]["q"], "1/64")

        self.assertEqual(
            [item["source_id"] for item in wb.pending_instructions()], ["D"]
        )
        self.assertEqual(wb.pending_instructions()[0]["token"], d_token)

        # 重开后幂等重试原登记：回执和原发令证据必须完全一致。
        wb.close()
        wb = open_workbook(self.path)
        retry = wb.record_execution(*submission)
        self.assertEqual(retry, receipt)
        self.assertEqual(wb.execution_count, 1)
        self.assertEqual(wb.remaining_capacity, 1)

        # 旧 B 指令已经撤销，不能登记。
        with self.assertRaises((ValueError, FollowUpConflict)):
            wb.record_execution(b_token, "shot-B-late")

        # 按验收剧本再次重开，然后恢复 A。
        wb.close()
        wb = open_workbook(self.path)

        # 恢复 A：B 优先于 D，因此 B 拿新令牌重新进入，D 被撤销。
        restored = wb.update_quality_flags(
            [
                {"id": "a0", "usable": True},
                {"id": "a1", "usable": True},
            ]
        )
        self.assertEqual(
            [item["source_id"] for item in restored["revoked"]], ["D"]
        )
        self.assertEqual(
            [item["source_id"] for item in restored["issued"]], ["B"]
        )
        new_b_token = restored["issued"][0]["token"]
        self.assertNotEqual(new_b_token, b_token)
        self.assertEqual(
            [item["source_id"] for item in wb.pending_instructions()], ["B"]
        )
        self.assertEqual(wb.pending_instructions()["B"]["token"], new_b_token)

        # 旧 B 令牌永久不可再用。
        with self.assertRaises((ValueError, FollowUpConflict)):
            wb.record_execution(b_token, "shot-B-old-token")

        # X 随后证据不足：合格总数 5 -> 4，B 的 q 再次变为 1/24 并离开计划，
        # D 顶替；但 X 的已执行记录不退回，容量余量仍只有一个。
        x_bad = [
            {"measurement_id": f"x{i}", "usable": False} for i in range(4)
        ]
        x_change = wb.update_quality_flags(x_bad)
        x_row = wb.candidate_results()["X"]
        self.assertEqual(x_row["status"], "insufficient_evidence")
        self.assertIsNone(x_row["q"])
        self.assertEqual(
            [item["source_id"] for item in x_change["revoked"]], ["B"]
        )
        self.assertEqual(x_change["revoked"][0]["token"], new_b_token)
        self.assertEqual(
            [item["source_id"] for item in x_change["issued"]], ["D"]
        )
        final_d_token = x_change["issued"][0]["token"]
        self.assertNotEqual(final_d_token, d_token)
        self.assertEqual(len(wb.executed_records()), 1)
        self.assertEqual(wb.executed_records()[0], receipt)
        self.assertEqual(wb.remaining_capacity, 1)
        self.assertEqual(len(wb.pending_instructions()), 1)
        self.assertEqual(
            wb.pending_instructions()[0]["source_id"], "D"
        )

        # 有效更正 + 未知 ID：整批拒绝，重开后一切如旧。
        before_tokens = [
            item["token"] for item in wb.pending_instructions()
        ]
        before_candidates = [dict(row) for row in wb.candidate_results()]
        bad_batch = [
            {"measurement_id": "b0", "usable": False},
            {"measurement_id": "does-not-exist", "usable": False},
        ]
        with self.assertRaises(ValueError):
            wb.update_quality_flags(bad_batch)

        wb.close()
        wb = open_workbook(self.path)
        self.assertEqual(
            [dict(row) for row in wb.candidate_results()],
            before_candidates,
        )
        b_row = wb.candidate_results()["B"]
        self.assertEqual((b_row["n"], b_row["q"]), (5, "1/24"))
        self.assertFalse(b_row["candidate"])
        after_tokens = [item["token"] for item in wb.pending_instructions()]
        self.assertEqual(after_tokens, before_tokens)
        self.assertEqual(
            [item["token"] for item in wb.pending_instructions()],
            [final_d_token],
        )
        self.assertEqual(wb.executed_records()[0], receipt)
        self.assertEqual(wb.execution_count, 1)
        self.assertEqual(wb.remaining_capacity, 1)

    def test_invalid_initialization_leaves_no_workbook(self):
        bad_path = Path(self.temp.name) / "missing" / "wb.json"
        batch = acceptance_batch()
        batch["fdr_threshold"] = "0"
        with self.assertRaises(ValueError):
            create_workbook(bad_path, batch)
        self.assertFalse(bad_path.exists())

    def test_duplicate_correction_pointer_rejected(self):
        with self.assertRaises(ValueError):
            self.wb.update_quality_flags(
                [
                    {"measurement_id": "a0", "usable": False},
                    {"measurement_id": "a0", "usable": True},
                ]
            )
        self.assertTrue(self.wb.candidate_results()["A"]["candidate"])


class ExactRationalAndOrdering(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def test_spelled_time_equality_and_inverse_variance_baseline_tie(self):
        path = Path(self.temp.name) / "rational.json"
        batch = {
            "sources": [
                {"id": "S", "baseline_magnitude": "10.4", "priority": 1}
            ],
            "measurements": [
                # 1/sigma^2 权重为 1 和 1/4，均值 = (10 + 12/4)/(1.25) = 10.4。
                {"id": "same-time-1", "source_id": "S", "mjd": "1.0",
                 "magnitude": "10", "sigma": "1", "usable": True},
                {"id": "same-time-2", "source_id": "S", "mjd": "1.000000",
                 "magnitude": "12", "sigma": "2.000000", "usable": True},
                {"id": "bright-1", "source_id": "S", "mjd": "2",
                 "magnitude": "9", "sigma": "0.1", "usable": True},
                {"id": "bright-2", "source_id": "S", "mjd": "3",
                 "magnitude": "9", "sigma": "0.1", "usable": True},
                {"id": "faint-1", "source_id": "S", "mjd": "4",
                 "magnitude": "11", "sigma": "0.1", "usable": True},
            ],
            "min_valid_times": 3,
            "fdr_threshold": "1",
            "capacity": 1,
        }
        wb = create_workbook(path, batch)
        row = wb.candidate_results()["S"]
        self.assertEqual(row["n"], 3)
        self.assertEqual(row["k"], 2)
        self.assertEqual(row["p"], "1/2")
        self.assertEqual(row["q"], "1/2")
        self.assertTrue(row["candidate"])

    def test_reordered_rows_do_not_change_results_or_tokens(self):
        first_path = Path(self.temp.name) / "first.json"
        second_path = Path(self.temp.name) / "second.json"
        batch1 = acceptance_batch()
        batch2 = copy.deepcopy(batch1)
        batch2["sources"].reverse()
        batch2["measurements"].reverse()

        first = create_workbook(first_path, batch1)
        second = create_workbook(second_path, batch2)
        self.assertEqual(
            first.candidate_results(), second.candidate_results()
        )
        self.assertEqual(
            [i["token"] for i in first.pending_instructions()],
            [i["token"] for i in second.pending_instructions()],
        )

    def test_retained_instruction_keeps_issue_evidence(self):
        path = Path(self.temp.name) / "retained.json"
        batch = acceptance_batch()
        # 给 B 增加一个偏亮时刻，使首次发令证据为 n=6；随后屏蔽一条变为 n=5，
        # 但 B 仍在计划内，原指令证据不得被当前统计量替换。
        batch["measurements"].append(
            {
                "id": "b5",
                "source_id": "B",
                "mjd": "205",
                "magnitude": "14.9",
                "sigma": "0.01",
                "usable": True,
            }
        )
        wb = create_workbook(path, batch)
        old_b = wb.pending_instructions()["B"]
        self.assertEqual(old_b["n"], 6)

        wb.update_quality_flags([{"id": "b0", "usable": False}])
        current_b = wb.candidate_results()["B"]
        self.assertEqual(current_b["n"], 5)
        retained_b = wb.pending_instructions()["B"]
        self.assertEqual(retained_b["token"], old_b["token"])
        self.assertEqual(retained_b["n"], 6)
        self.assertEqual(retained_b["p"], old_b["p"])
        self.assertEqual(retained_b["q"], old_b["q"])
        self.assertNotEqual(retained_b["p"], current_b["p"])

    def test_p_ties_share_q_and_q_equal_threshold_is_candidate(self):
        path = Path(self.temp.name) / "ties.json"
        sources = [
            {"id": "b", "baseline_magnitude": "15", "priority": 1},
            {"id": "a", "baseline_magnitude": "15", "priority": 2},
        ]
        measurements = []
        for sid in ("a", "b"):
            for index in range(5):
                measurements.append(
                    {
                        "id": f"{sid}{index}",
                        "source_id": sid,
                        "mjd": str(index),
                        "magnitude": "14",
                        "sigma": "0.1",
                        "usable": True,
                    }
                )
        wb = create_workbook(
            path,
            {
                "sources": sources,
                "measurements": measurements,
                "min_valid_times": 5,
                # 1/32 = 0.03125，恰好等于并列 BH q。
                "fdr_threshold": "0.03125",
                "capacity": 1,
            },
        )
        rows = wb.candidate_results()
        self.assertEqual(rows["a"]["p"], "1/32")
        self.assertEqual(rows["b"]["p"], "1/32")
        self.assertEqual(rows["a"]["q"], "1/32")
        self.assertEqual(rows["b"]["q"], "1/32")
        self.assertTrue(rows["a"]["candidate"])
        self.assertTrue(rows["b"]["candidate"])
        self.assertEqual(
            [i["source_id"] for i in wb.pending_instructions()], ["b"]
        )


class RegistrationConflicts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "conflicts.json"
        self.wb = create_workbook(self.path, acceptance_batch())

    def test_unknown_and_cross_identity_conflicts_persist(self):
        pending = self.wb.pending_instructions()
        x_token = pending["X"]["token"]
        b_token = pending["B"]["token"]
        receipt = self.wb.record_execution(x_token, "shot-1")

        with self.assertRaises(FollowUpConflict):
            self.wb.record_execution(x_token, "shot-2")
        with self.assertRaises(FollowUpConflict):
            self.wb.record_execution(b_token, "shot-1")
        with self.assertRaises(FollowUpConflict):
            self.wb.record_execution("fu-999999999999", "shot-3")

        # 失败登记不改变待执行列表和容量。
        self.assertEqual(
            [i["source_id"] for i in self.wb.pending_instructions()], ["B"]
        )
        self.assertEqual(self.wb.execution_count, 1)

        reopened = open_workbook(self.path)
        self.assertEqual(reopened.executed_records()[0], receipt)
        with self.assertRaises(FollowUpConflict):
            reopened.record_execution(b_token, "shot-1")


class InputBounds(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def test_source_without_measurements_is_allowed(self):
        path = Path(self.temp.name) / "empty-source.json"
        wb = create_workbook(
            path,
            {
                "sources": [
                    {"id": "EMPTY", "baseline_magnitude": "15", "priority": 1},
                    {"id": "USED", "baseline_magnitude": "15", "priority": 2},
                ],
                "measurements": [
                    {"id": "u1", "source_id": "USED", "mjd": "1",
                     "magnitude": "14", "sigma": "0.1", "usable": True}
                ],
                "min_valid_times": 3,
                "fdr_threshold": "0.04",
                "capacity": 1,
            },
        )
        row = wb.candidate_results()["EMPTY"]
        self.assertEqual(row["n"], 0)
        self.assertEqual(row["k"], 0)
        self.assertEqual(row["status"], "insufficient_evidence")

    def test_max_128_distinct_moments(self):
        def batch(count):
            return {
                "sources": [
                    {"id": "S", "baseline_magnitude": "15", "priority": 1}
                ],
                "measurements": [
                    {"id": f"m{i}", "source_id": "S", "mjd": str(i),
                     "magnitude": "14.9", "sigma": "0.1", "usable": True}
                    for i in range(count)
                ],
                "min_valid_times": 3,
                "fdr_threshold": "1",
                "capacity": 1,
            }

        ok_path = Path(self.temp.name) / "ok.json"
        wb = create_workbook(ok_path, batch(128))
        self.assertEqual(wb.candidate_results()["S"]["n"], 128)

        bad_path = Path(self.temp.name) / "bad.json"
        with self.assertRaises(ValueError):
            create_workbook(bad_path, batch(129))
        self.assertFalse(bad_path.exists())


if __name__ == "__main__":
    unittest.main()
