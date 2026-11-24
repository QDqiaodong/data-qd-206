import copy
import json
import math
import random
import unittest

from star_crossmatch import solve_astrometry


def tangent_basis(ra0=0.0, dec0=0.0):
    r0 = math.radians(ra0)
    d0 = math.radians(dec0)
    return r0, d0


def inverse_gnomonic(xi_deg, eta_deg, ra0=0.0, dec0=0.0):
    """测试侧独立实现的标准 gnomonic 反投影。"""
    ra0r = math.radians(ra0)
    dec0r = math.radians(dec0)
    xi = math.radians(xi_deg)
    eta = math.radians(eta_deg)
    rho = math.hypot(xi, eta)
    if rho == 0.0:
        return ra0, dec0
    c = math.atan(rho)
    dec = math.asin(
        math.cos(c) * math.sin(dec0)
        + eta * math.sin(c) * math.cos(dec0) / rho
    )
    ra = ra0r + math.atan2(
        xi * math.sin(c),
        rho * math.cos(dec0r) * math.cos(c)
        - eta * math.sin(dec0r) * math.sin(c),
    )
    return math.degrees(ra) % 360.0, math.degrees(dec)


def sky_for(x, y, affine, ra0=0.0, dec0=0.0):
    a, b, c, d, e, f = affine
    return inverse_gnomonic(a * x + b * y + c, d * x + e * y + f, ra0, dec0)


def decimal6(value):
    return f"{value:.6f}"


def control(cid, x, y, ra, dec):
    return {
        "id": cid,
        "x": str(x) if isinstance(x, str) else decimal6(x),
        "y": str(y) if isinstance(y, str) else decimal6(y),
        "ra": decimal6(ra),
        "dec": decimal6(dec),
    }


def detection(did, x, y):
    return {"id": did, "x": decimal6(x), "y": decimal6(y)}


def valid_data():
    affine = (0.010, 0.002, 0.050, 0.001, 0.012, -0.030)
    grid = [(x, y) for x in range(4) for y in range(4)]
    controls = []
    for i, (x, y) in enumerate(grid):
        ra, dec = sky_for(x, y, affine)
        # 约 0.03 角秒的确定性小噪声，仍远小于 0.5 角秒上限。
        noise_ra = 0.000002 * ((i % 3) - 1)
        noise_dec = 0.000002 * ((i % 2) * 2 - 1)
        controls.append(control(
            f"g{i:02d}", x, y, ra + noise_ra, dec + noise_dec
        ))
    controls.append(control("bad1", 7.5, 7.5, 12.0, 8.0))
    controls.append(control("bad2", -4.0, 6.0, 350.0, -5.0))
    detections = [
        detection("inside", 1.5, 1.5),
        detection("edge", 3.0, 1.5),
        detection("outside", 12.0, 12.0),
    ]
    return {
        "tangent_point": {"ra": "0", "dec": "0"},
        "controls": controls,
        "detections": detections,
        "residual_limit_arcsec": "0.5",
        "min_inliers": 6,
    }, affine


class SolvedAstrometry(unittest.TestCase):
    def test_known_affine_noise_and_bad_controls(self):
        data, affine = valid_data()
        original = copy.deepcopy(data)
        result = solve_astrometry(data)
        self.assertEqual(data, original)  # 不修改输入

        self.assertEqual(result["status"], "solved")
        self.assertTrue(result["valid"])
        inliers = {c["id"]: c["inlier"] for c in result["controls"]}
        self.assertEqual(len([v for v in inliers.values() if v]), 16)
        self.assertFalse(inliers["bad1"])
        self.assertFalse(inliers["bad2"])

        params = result["affine"]
        self.assertAlmostEqual(float(params["a"]), affine[0], delta=2e-6)
        self.assertAlmostEqual(float(params["b"]), affine[1], delta=2e-6)
        self.assertAlmostEqual(float(params["c"]), affine[2], delta=2e-6)
        self.assertAlmostEqual(float(params["d"]), affine[3], delta=2e-6)
        self.assertAlmostEqual(float(params["e"]), affine[4], delta=2e-6)
        self.assertAlmostEqual(float(params["f"]), affine[5], delta=2e-6)

        rows = {d["id"]: d for d in result["detections"]}
        self.assertEqual(rows["outside"]["status"], "outside_footprint")
        self.assertIsNone(rows["outside"]["ra"])
        self.assertIsNone(rows["outside"]["dec"])
        for did in ("inside", "edge"):
            self.assertIsNotNone(rows[did]["ra"])
            x = float(data["detections"][[
                d["id"] for d in data["detections"]
            ].index(did)]["x"])
            y = float(data["detections"][[
                d["id"] for d in data["detections"]
            ].index(did)]["y"])
            ra, dec = sky_for(x, y, affine)
            self.assertAlmostEqual(float(rows[did]["ra"]), ra, delta=0.002)
            self.assertAlmostEqual(float(rows[did]["dec"]), dec, delta=0.002)

        # JSON 可序列化且不含 NaN/Infinity。
        json.dumps(result, allow_nan=False, sort_keys=True)

    def test_boundary_vertex_removal_moves_detection_outside_hull(self):
        affine = (0.01, 0.0, 0.1, 0.0, 0.01, 0.1)
        square = [(0, 0), (10, 0), (10, 10), (0, 10)]
        controls = []
        for cid, (x, y) in zip(("bl", "br", "tr", "tl"), square):
            ra, dec = sky_for(x, y, affine)
            controls.append(control(cid, x, y, ra, dec))
        data = {
            "tangent_point": {"ra": "0", "dec": "0"},
            "controls": controls,
            "detections": [detection("edge-mid", 10, 5)],
            "residual_limit_arcsec": "0.05",
            "min_inliers": 3,
        }
        first = solve_astrometry(data)
        self.assertEqual(first["status"], "solved")
        det = first["detections"][0]
        self.assertTrue(det["inside_footprint"])
        self.assertEqual(det["status"], "transformed")
        self.assertIsNotNone(det["ra"])

        # 删除决定右上角/凸包边界的控制星后，剩余三角形为 x+y<=10；
        # (10,5) 明确落在凸包外，不能外推。
        data["controls"] = [c for c in controls if c["id"] != "tr"]
        second = solve_astrometry(data)
        det = second["detections"][0]
        self.assertFalse(det["inside_footprint"])
        self.assertEqual(det["status"], "outside_footprint")
        self.assertIsNone(det["ra"])
        self.assertIsNone(det["dec"])

    def test_ra_zero_wrap_projection_and_inverse(self):
        affine = (0.01, 0.0, 0.0, 0.0, 0.01, 0.0)
        positions = [(-5, -1), (5, -1), (-5, 1), (5, 1)]
        controls = []
        for i, (x, y) in enumerate(positions):
            ra, dec = sky_for(x, y, affine)
            controls.append(control(f"wrap{i}", x, y, ra, dec))
        data = {
            "tangent_point": {"ra": "0", "dec": "0"},
            "controls": controls,
            "detections": [detection("cross", -2.5, 0)],
            "residual_limit_arcsec": "0.05",
            "min_inliers": 3,
        }
        result = solve_astrometry(data)
        self.assertEqual(result["status"], "solved")
        row = result["detections"][0]
        ra, dec = sky_for(-2.5, 0, affine)
        self.assertGreater(float(row["ra"]), 359.0)
        self.assertLess(float(row["ra"]), 360.0)
        self.assertAlmostEqual(float(row["ra"]), ra, delta=1e-7)
        self.assertAlmostEqual(float(row["dec"]), dec, delta=1e-8)
        projected = {c["id"]: c for c in result["controls"]}
        self.assertLess(float(projected["wrap0"]["xi"]), 0.0)


class FailureAndTie(unittest.TestCase):
    def base(self, controls, min_inliers=3, limit="0.05"):
        return {
            "tangent_point": {"ra": "0", "dec": "0"},
            "controls": controls,
            "detections": [detection("d", 0.5, 0.5)],
            "residual_limit_arcsec": limit,
            "min_inliers": min_inliers,
        }

    def test_collinear_controls_rank_deficient(self):
        controls = []
        for i in range(5):
            ra, dec = sky_for(i, 2 * i, (0.01, 0.02, 0.0, 0.005, 0.01, 0.0))
            controls.append(control(f"line{i}", i, 2 * i, ra, dec))
        result = solve_astrometry(self.base(controls))
        self.assertEqual(result["status"], "rank_deficient")
        self.assertIsNone(result["affine"])
        self.assertFalse(result["valid"])
        self.assertEqual(result["detections"][0]["status"], "rank_deficient")
        self.assertIsNone(result["detections"][0]["ra"])
        self.assertTrue(all(c["xi"] is not None for c in result["controls"]))

    def test_insufficient_inliers_reported_separately(self):
        controls = [
            control("a", 0, 0, 0.0, 0.0),
            control("b", 1, 0, 0.01, 0.0),
            control("c", 0, 1, 0.0, 0.01),
            control("d", 1, 1, 4.0, 4.0),
            control("e", 2, 2, 357.0, 2.0),
        ]
        result = solve_astrometry(self.base(controls, min_inliers=4))
        self.assertEqual(result["status"], "insufficient_inliers")
        self.assertIsNone(result["affine"])
        self.assertEqual(result["inlier_count"], 3)
        self.assertTrue(all(d["ra"] is None for d in result["detections"]))

    def test_two_equal_best_consensus_solutions_are_ambiguous(self):
        shape = [(0, 0), (1, 0), (0, 1), (1, 1)]
        affine_a = (0.02, 0.0, 0.0, 0.0, 0.02, 0.0)
        affine_b = (-0.02, 0.0, 0.0, 0.0, 0.02, 0.0)
        controls = []
        for prefix, affine in (("A", affine_a), ("B", affine_b)):
            for i, (x, y) in enumerate(shape):
                ra, dec = sky_for(x, y, affine)
                controls.append(control(f"{prefix}{i}", x, y, ra, dec))
        result = solve_astrometry(self.base(controls, min_inliers=4, limit="0.1"))
        self.assertEqual(result["status"], "ambiguous_solution")
        self.assertIsNone(result["affine"])
        self.assertGreaterEqual(result["tied_solution_count"], 2)
        self.assertGreater(
            float(result["maximum_detection_difference_arcsec"]), 0.1
        )
        row = result["detections"][0]
        self.assertEqual(row["status"], "ambiguous_solution")
        self.assertIsNone(row["ra"])
        self.assertIsNone(row["dec"])
        self.assertGreaterEqual(len(result["alternatives"]), 2)

    def test_reorder_rows_preserves_identity_based_json_result(self):
        data, _ = valid_data()
        r1 = solve_astrometry(data)
        data2 = copy.deepcopy(data)
        random.Random(206).shuffle(data2["controls"])
        random.Random(207).shuffle(data2["detections"])
        r2 = solve_astrometry(data2)
        j1 = json.dumps(r1, sort_keys=True, ensure_ascii=False)
        j2 = json.dumps(r2, sort_keys=True, ensure_ascii=False)
        self.assertEqual(json.loads(j2), json.loads(j1))
        self.assertEqual([c["id"] for c in r2["controls"]],
                         sorted(c["id"] for c in data["controls"]))
        self.assertEqual([d["id"] for d in r2["detections"]],
                         sorted(d["id"] for d in data["detections"]))


class InvalidInput(unittest.TestCase):
    def setUp(self):
        self.data, _ = valid_data()

    def expect_value_error(self, mutate):
        data = copy.deepcopy(self.data)
        mutate(data)
        before = copy.deepcopy(data)
        with self.assertRaises(ValueError):
            solve_astrometry(data)
        self.assertEqual(data, before)

    def test_strict_fields_and_types(self):
        self.expect_value_error(lambda d: d.pop("min_inliers"))
        self.expect_value_error(lambda d: d.update(extra="x"))
        self.expect_value_error(lambda d: d["tangent_point"].update(extra=1))
        self.expect_value_error(lambda d: d["controls"][0].pop("ra"))
        self.expect_value_error(lambda d: d["controls"][0].update(extra=1))
        self.expect_value_error(lambda d: d["detections"][0].pop("x"))
        self.expect_value_error(lambda d: d.update(min_inliers="6"))
        self.expect_value_error(lambda d: d.update(min_inliers=True))
        self.expect_value_error(lambda d: d.update(residual_limit_arcsec=0.5))
        self.expect_value_error(lambda d: d.update(residual_limit_arcsec="0"))

    def test_decimal_syntax_and_ranges(self):
        self.expect_value_error(lambda d: d["controls"][0].update(x=1))
        self.expect_value_error(lambda d: d["controls"][0].update(x=".5"))
        self.expect_value_error(lambda d: d["controls"][0].update(y="1e3"))
        self.expect_value_error(lambda d: d["controls"][0].update(ra="360"))
        self.expect_value_error(lambda d: d["controls"][0].update(dec="91"))
        self.expect_value_error(
            lambda d: d["tangent_point"].update(ra="360.000001")
        )

    def test_identities_duplicates_and_bounds(self):
        self.expect_value_error(
            lambda d: d["controls"][0].update(id="bad id")
        )
        self.expect_value_error(
            lambda d: d["controls"][1].update(id=d["controls"][0]["id"])
        )
        self.expect_value_error(
            lambda d: d["detections"][1].update(id=d["detections"][0]["id"])
        )
        self.expect_value_error(lambda d: d.update(min_inliers=2))
        self.expect_value_error(lambda d: d.update(min_inliers=1001))
        self.expect_value_error(lambda d: d.update(controls=[]))
        self.expect_value_error(lambda d: d.update(detections=[]))

    def test_hemisphere_rejected(self):
        def move(d):
            d["controls"][0]["ra"] = "90"
            d["controls"][0]["dec"] = "0"
            d["tangent_point"] = {"ra": "0", "dec": "0"}
        self.expect_value_error(move)


if __name__ == "__main__":
    unittest.main()
