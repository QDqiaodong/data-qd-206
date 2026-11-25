"""``infer_distance_catalog`` 的离线验收测试。

覆盖：

* 固定小批次对照独立高精度 oracle（保存在 distance_oracle_fixture.json，
  由 tests/_distance_oracle.py 这套独立 Gauss-Legendre 实现离线生成）；
* 负/零/正视差、退化零点先验、极小/极大标准差与尺度、严格为正且有序的分位数；
* 共享零点是全批耦合的公共随机量（不是逐源独立误差），第三源会改变前两源；
* pair 反向概率互补、多 pair 共用源、输入行置换不改变按身份比较的结果；
* 严格字段/类型/身份/范围校验，失败整体抛 ValueError 且不改输入。

所有测试只用标准库、不访问网络。
"""
import copy
import json
import math
import os
import unittest

from star_crossmatch import (
    DISTANCE_CATALOG_FORMAT,
    infer_distance_catalog,
)

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_PATH = os.path.join(HERE, "distance_oracle_fixture.json")

DIST_TOL = 1e-4
PAIR_TOL = 1e-5


def fixture_input():
    with open(FIXTURE_PATH, encoding="utf-8") as handle:
        return json.load(handle)


class OracleConsistency(unittest.TestCase):
    """与独立确定性积分 oracle 的固定数值对照。"""

    @classmethod
    def setUpClass(cls):
        with open(FIXTURE_PATH, encoding="utf-8") as handle:
            cls.bundle = json.load(handle)
        cls.oracle = cls.bundle["oracle"]
        cls.result = infer_distance_catalog(copy.deepcopy(cls.bundle["input"]))

    def test_format_version_and_finite_json(self):
        self.assertEqual(self.result["format"], DISTANCE_CATALOG_FORMAT)
        # 全部为有限 JSON 数字（allow_nan=False 会拒绝 NaN/Infinity）。
        json.dumps(self.result, allow_nan=False)

    def test_zero_point_summary_matches_oracle(self):
        zp = self.result["zero_point"]
        ozp = self.oracle["zero_point"]
        self.assertAlmostEqual(
            zp["mean_mas"], ozp["mean_mas"], delta=DIST_TOL)
        # 零点均值/标准差本身是 mas，量级远小于 pc，用更紧的绝对容差核对。
        self.assertAlmostEqual(
            zp["sigma_mas"], ozp["sigma_mas"], delta=1e-6)

    def test_source_quantiles_match_oracle(self):
        rows = {row["id"]: row for row in self.result["sources"]}
        for sid, expected in self.oracle["sources"].items():
            got = rows[sid]
            self.assertAlmostEqual(
                got["median_pc"], expected["median_pc"], delta=DIST_TOL)
            self.assertAlmostEqual(
                got["ci_lower_pc"], expected["ci_lower_pc"],
                delta=DIST_TOL)
            self.assertAlmostEqual(
                got["ci_upper_pc"], expected["ci_upper_pc"],
                delta=DIST_TOL)

    def test_pair_probability_matches_oracle(self):
        pairs = {row["id"]: row for row in self.result["pairs"]}
        expected = self.oracle["pairs"]["A-before-B"]
        self.assertAlmostEqual(
            pairs["A-before-B"]["probability"], expected, delta=PAIR_TOL)

    def test_oracle_engine_is_independent_of_production(self):
        """oracle 模块不得导入被测包，保证是真正独立的实现。"""
        import _distance_oracle as oracle  # noqa: PLC0415
        source_path = oracle.__file__
        with open(source_path, encoding="utf-8") as handle:
            text = handle.read()
        self.assertNotIn("infer_distance_catalog", text)
        self.assertNotIn("star_crossmatch", text)
        # 同一输入下两套实现的距离分位数确实分别落在容差内（交叉验证）。
        independent = oracle.run(copy.deepcopy(self.bundle["input"]))
        rows = {row["id"]: row for row in self.result["sources"]}
        for sid in ("A", "B"):
            self.assertAlmostEqual(
                rows[sid]["median_pc"],
                independent["sources"][sid]["median_pc"],
                delta=DIST_TOL)


class PositiveOrderedQuantiles(unittest.TestCase):
    """负、零、正视差都给出有限、严格为正且有序的距离分位数。"""

    def _run(self, parallax, sigma, zp, L="1200"):
        data = {
            "scale_length_pc": L,
            "zero_point": zp,
            "sources": [{
                "id": "S",
                "parallax_mas": parallax,
                "sigma_mas": sigma,
            }],
            "pairs": [],
        }
        return infer_distance_catalog(data)

    def _assert_well_formed(self, result):
        row = result["sources"][0]
        for key in ("median_pc", "ci_lower_pc", "ci_upper_pc"):
            self.assertTrue(math.isfinite(row[key]), key)
            self.assertGreater(row[key], 0.0, key)
        self.assertLess(row["ci_lower_pc"], row["median_pc"])
        self.assertLess(row["median_pc"], row["ci_upper_pc"])

    def test_negative_parallax_fixed_zero(self):
        result = self._run("-0.20", "0.50",
                           {"mean_mas": "0.08", "sigma_mas": "0"})
        self._assert_well_formed(result)
        self.assertEqual(result["zero_point"]["sigma_mas"], 0.0)
        self.assertEqual(result["zero_point"]["mean_mas"], 0.08)

    def test_negative_parallax_uncertain_zero(self):
        # 负视差 + 共享零点不确定性：后验完全由先验与高斯似然共同决定，
        # 仍给出有限、严格为正且有序的距离。
        result = self._run("-0.80", "0.25",
                           {"mean_mas": "0.0", "sigma_mas": "0.2"})
        self._assert_well_formed(result)
        self.assertGreater(result["zero_point"]["sigma_mas"], 0.0)

    def test_zero_parallax(self):
        result = self._run("0", "0.01",
                           {"mean_mas": "0", "sigma_mas": "0.1"})
        self._assert_well_formed(result)

    def test_large_positive_parallax_not_naive_inverse(self):
        result = self._run("50", "0.01",
                           {"mean_mas": "0", "sigma_mas": "0"})
        self._assert_well_formed(result)
        median = result["sources"][0]["median_pc"]
        # 朴素倒数给 20 pc；贝叶斯后验中位数接近但不应恰好等于倒数。
        self.assertAlmostEqual(median, 20.0, delta=1.0)
        self.assertNotEqual(median, 1000.0 / 50.0)
        # 区间端点也不得取倒数。
        self.assertNotEqual(result["sources"][0]["ci_lower_pc"], 20.0)
        self.assertNotEqual(result["sources"][0]["ci_upper_pc"], 20.0)

    def test_very_large_and_very_small_sigma(self):
        big = self._run("1.0", "9999",
                        {"mean_mas": "0", "sigma_mas": "0"})
        self._assert_well_formed(big)
        small = self._run("0.5", "0.00001",
                          {"mean_mas": "0", "sigma_mas": "0"})
        self._assert_well_formed(small)

    def test_extreme_scale_lengths(self):
        near = self._run("5.0", "0.01",
                         {"mean_mas": "0", "sigma_mas": "0"}, L="0.01")
        far = self._run("0.001", "0.0005",
                        {"mean_mas": "0", "sigma_mas": "0"}, L="1000000")
        self._assert_well_formed(near)
        self._assert_well_formed(far)

    def test_negative_zero_point_mean(self):
        result = self._run("-2.0", "0.3",
                           {"mean_mas": "-1.5", "sigma_mas": "0.2"})
        self._assert_well_formed(result)


class SharedZeroPointCoupling(unittest.TestCase):
    """共享零点是公共随机量，所有源共同更新其后验。"""

    def _two_similar_sources(self, informed):
        sources = [
            {"id": "S1", "parallax_mas": "0.40", "sigma_mas": "0.30"},
            {"id": "S2", "parallax_mas": "0.42", "sigma_mas": "0.28"},
        ]
        if informed:
            # 第三源测量精确，对零点有强信息量。
            sources.append(
                {"id": "S3", "parallax_mas": "1.20", "sigma_mas": "0.01"})
        return {
            "scale_length_pc": "1200",
            "zero_point": {"mean_mas": "0.00", "sigma_mas": "0.20"},
            "sources": sources,
            "pairs": [{"id": "p12", "nearer": "S1", "farther": "S2"}],
        }

    def test_informative_third_source_moves_zero_and_distances(self):
        before = infer_distance_catalog(self._two_similar_sources(False))
        after = infer_distance_catalog(self._two_similar_sources(True))
        # 精确第三源参与联合更新，改变共享零点的后验均值与标准差。
        self.assertNotAlmostEqual(
            before["zero_point"]["mean_mas"],
            after["zero_point"]["mean_mas"], places=6)
        self.assertNotAlmostEqual(
            before["zero_point"]["sigma_mas"],
            after["zero_point"]["sigma_mas"], places=6)
        # 前两源（不把第三源当定标标准，只通过共享 z 耦合）的距离随之改变。
        b = {r["id"]: r for r in before["sources"]}
        a = {r["id"]: r for r in after["sources"]}
        moved = any(
            abs(b[sid]["median_pc"] - a[sid]["median_pc"]) > 1e-3
            for sid in ("S1", "S2"))
        self.assertTrue(moved)

    def test_joint_differs_from_fixed_prior_mean(self):
        # 把零点固定在先验均值（sigma=0）应得到与联合后验不同的结果。
        joint = infer_distance_catalog(self._two_similar_sources(True))
        fixed_data = self._two_similar_sources(True)
        fixed_data["zero_point"] = {"mean_mas": "0.00", "sigma_mas": "0"}
        fixed = infer_distance_catalog(fixed_data)
        self.assertNotAlmostEqual(
            joint["zero_point"]["mean_mas"],
            fixed["zero_point"]["mean_mas"], places=6)

    def test_joint_differs_from_independent_per_source_zero(self):
        """逐源独立抽取零点等价于把零点不确定性并入各源、源间解耦。

        用 sigma=0 的固定零点逐源单独推断，再与全批联合结果比较：联合模型
        因共享 z 而在源间传递信息，二者的零点边缘与距离不应相同。
        """
        data = self._two_similar_sources(True)
        joint = infer_distance_catalog(data)
        independent_medians = {}
        for src in data["sources"]:
            single = {
                "scale_length_pc": data["scale_length_pc"],
                # 每源各自面对一个固定在先验均值的零点，源间无耦合。
                "zero_point": {"mean_mas": "0.00", "sigma_mas": "0"},
                "sources": [copy.deepcopy(src)],
                "pairs": [],
            }
            out = infer_distance_catalog(single)
            independent_medians[src["id"]] = out["sources"][0]["median_pc"]
        joint_medians = {r["id"]: r["median_pc"]
                         for r in joint["sources"]}
        self.assertTrue(any(
            abs(joint_medians[sid] - independent_medians[sid]) > 1e-3
            for sid in joint_medians))


class PairOrdering(unittest.TestCase):
    def setUp(self):
        self.data = fixture_input()["input"]

    def test_reversed_pair_probabilities_sum_to_one(self):
        forward = infer_distance_catalog(copy.deepcopy(self.data))
        reversed_data = copy.deepcopy(self.data)
        reversed_data["pairs"] = [{
            "id": "rev", "nearer": "B", "farther": "A"}]
        backward = infer_distance_catalog(reversed_data)
        p_forward = forward["pairs"][0]["probability"]
        p_backward = backward["pairs"][0]["probability"]
        self.assertTrue(0.0 <= p_forward <= 1.0)
        self.assertTrue(0.0 <= p_backward <= 1.0)
        self.assertAlmostEqual(p_forward + p_backward, 1.0,
                               delta=PAIR_TOL)

    def test_multiple_pairs_share_sources(self):
        data = copy.deepcopy(self.data)
        data["pairs"] = [
            {"id": "pAB", "nearer": "A", "farther": "B"},
            {"id": "pBA", "nearer": "B", "farther": "A"},
            {"id": "pAA-x", "nearer": "A", "farther": "B"},
        ]
        result = infer_distance_catalog(data)
        self.assertEqual(len(result["pairs"]), 3)
        ids = [p["id"] for p in result["pairs"]]
        self.assertEqual(ids, sorted(ids))  # 按 pair ID 排序
        probs = {p["id"]: p["probability"] for p in result["pairs"]}
        self.assertAlmostEqual(
            probs["pAB"] + probs["pBA"], 1.0, delta=PAIR_TOL)

    def test_empty_pairs_allowed(self):
        data = copy.deepcopy(self.data)
        data["pairs"] = []
        result = infer_distance_catalog(data)
        self.assertEqual(result["pairs"], [])


class OrderIndependence(unittest.TestCase):
    def test_row_permutation_identity_equal_json(self):
        data = fixture_input()["input"]
        first = infer_distance_catalog(copy.deepcopy(data))
        permuted = copy.deepcopy(data)
        permuted["sources"] = list(reversed(permuted["sources"]))
        permuted["pairs"] = list(reversed(permuted["pairs"]))
        second = infer_distance_catalog(permuted)
        self.assertEqual(
            json.dumps(first, sort_keys=True),
            json.dumps(second, sort_keys=True))
        # 输出始终按身份排序。
        source_ids = [r["id"] for r in first["sources"]]
        self.assertEqual(source_ids, sorted(source_ids))
        pair_ids = [r["id"] for r in first["pairs"]]
        self.assertEqual(pair_ids, sorted(pair_ids))


class StrictValidation(unittest.TestCase):
    def setUp(self):
        self.data = fixture_input()["input"]

    def expect_value_error(self, mutate):
        data = copy.deepcopy(self.data)
        mutate(data)
        before = copy.deepcopy(data)
        with self.assertRaises(ValueError):
            infer_distance_catalog(data)
        # 失败不得修改调用者输入，也不返回部分结果。
        self.assertEqual(data, before)

    def test_top_level_fields(self):
        self.expect_value_error(lambda d: d.pop("sources"))
        self.expect_value_error(lambda d: d.pop("pairs"))
        self.expect_value_error(lambda d: d.update(extra=1))
        self.expect_value_error(lambda d: d.__setitem__("scale_length_pc", 1200))

    def test_nested_fields(self):
        self.expect_value_error(
            lambda d: d["zero_point"].update(extra=1))
        self.expect_value_error(
            lambda d: d["zero_point"].pop("sigma_mas"))
        self.expect_value_error(
            lambda d: d["sources"][0].pop("sigma_mas"))
        self.expect_value_error(
            lambda d: d["sources"][0].update(extra=1))
        self.expect_value_error(
            lambda d: d["pairs"][0].update(extra=1))

    def test_types_and_bool_impersonation(self):
        self.expect_value_error(
            lambda d: d.__setitem__("sources", "not-a-list"))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("parallax_mas", -0.2))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("parallax_mas", True))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("sigma_mas", False))
        self.expect_value_error(
            lambda d: d["zero_point"].__setitem__("mean_mas", 0.08))
        self.expect_value_error(
            lambda d: d["zero_point"].__setitem__("sigma_mas", True))

    def test_decimal_syntax(self):
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("parallax_mas", ".5"))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("parallax_mas", "1e3"))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("sigma_mas", "nan"))
        self.expect_value_error(
            lambda d: d.__setitem__("scale_length_pc", "abc"))

    def test_ranges(self):
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("sigma_mas", "0"))
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("sigma_mas", "-1"))
        self.expect_value_error(
            lambda d: d.__setitem__("scale_length_pc", "0"))
        self.expect_value_error(
            lambda d: d.__setitem__("scale_length_pc", "-1200"))
        self.expect_value_error(
            lambda d: d["zero_point"].__setitem__("sigma_mas", "-0.1"))

    def test_identities_and_endpoints(self):
        # 重复 source id
        self.expect_value_error(lambda d: d["sources"][1].__setitem__(
            "id", d["sources"][0]["id"]))
        # 非法 id 字符
        self.expect_value_error(
            lambda d: d["sources"][0].__setitem__("id", "bad id"))
        # 未知 pair 端点
        self.expect_value_error(
            lambda d: d["pairs"][0].__setitem__("farther", "GHOST"))
        # 两端相同
        self.expect_value_error(lambda d: d["pairs"][0].__setitem__(
            "farther", d["pairs"][0]["nearer"]))
        # 重复 pair id
        self.expect_value_error(lambda d: d["pairs"].append(
            {"id": d["pairs"][0]["id"], "nearer": "B", "farther": "A"}))

    def test_array_bounds(self):
        self.expect_value_error(lambda d: d.__setitem__("sources", []))
        too_many = [
            {"id": f"S{i}", "parallax_mas": "0.1", "sigma_mas": "0.2"}
            for i in range(65)]
        self.expect_value_error(lambda d: d.__setitem__("sources", too_many))
        too_many_pairs = [
            {"id": f"P{i}", "nearer": "A", "farther": "B"}
            for i in range(257)]
        self.expect_value_error(
            lambda d: d.__setitem__("pairs", too_many_pairs))


if __name__ == "__main__":
    unittest.main()
