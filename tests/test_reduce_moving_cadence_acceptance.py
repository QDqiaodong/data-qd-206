import copy
import json
import unittest

from star_crossmatch import process_exposures, reduce_moving_cadence


def obs(oid, ra="0", dec="0", mag="13", sigma="0.001"):
    return {
        "id": oid,
        "ra": ra,
        "dec": dec,
        "sigma": sigma,
        "instrumental_magnitude": mag,
    }


def ref(rid, ra, mag="15", calibrator=True):
    return {
        "id": rid,
        "ra": ra,
        "dec": "0",
        "magnitude": mag,
        "calibrator": calibrator,
    }


def position(eid, oid, x, y="0"):
    return {"exposure": eid, "observation": oid, "x": str(x), "y": str(y)}


def first_cadence_data(middle_mjd="60001"):
    photometry = {
        "reference": [ref("S1", "10"), ref("S2", "20"), ref("S3", "30")],
        "threshold": "1",
        "anchor_tolerance": "0.2",
        "exposures": [
            {"id": "E", "observations": [obs("A", ra="0")]},
            {
                "id": "M",
                "observations": [
                    obs("B", ra="30"),
                    obs("D", ra="40"),
                    obs("P1", ra="10"),
                    obs("P2", ra="20"),
                    obs("F", ra="50"),
                ],
            },
            {"id": "L", "observations": [obs("C", ra="1")]},
        ],
    }
    return {
        "photometry": photometry,
        "epochs": [
            {"exposure": "E", "mjd": "60000"},
            {"exposure": "M", "mjd": middle_mjd},
            {"exposure": "L", "mjd": "60002"},
        ],
        "positions": [
            position("E", "A", "0"),
            position("M", "B", "0.9"),
            position("M", "D", "1.1"),
            position("M", "P1", "10", "10"),
            position("M", "P2", "11", "10"),
            position("M", "F", "50", "20"),
            position("L", "C", "2"),
        ],
        "min_speed": "0.75",
        "max_speed": "1.25",
        "middle_tolerance": "0.2",
    }


def record(result, eid, oid):
    matches = [
        item
        for item in result["detections"]
        if item["exposure"] == eid and item["observation"] == oid
    ]
    assert len(matches) == 1
    return matches[0]


def identity_keys(triple):
    return [(item["exposure"], item["observation"]) for item in triple]


class MovingCadenceAcceptanceTest(unittest.TestCase):
    def test_standard_star_overlap_and_counterfactuals(self):
        data = first_cadence_data()
        original = copy.deepcopy(data)
        result = reduce_moving_cadence(data)
        self.assertEqual(data, original)
        self.assertEqual(json.loads(json.dumps(result)), result)

        self.assertEqual(len(result["components"]), 1)
        component = result["components"][0]
        self.assertEqual(component["detection_count"], 4)
        self.assertEqual(component["unconstrained_best"]["trajectory_count"], 1)
        self.assertEqual(component["unconstrained_best"]["total_cost"], "1/100")
        self.assertEqual(result["global_best"]["total_cost"], "1/100")

        self.assertEqual(record(result, "E", "A")["state"], "always")
        self.assertEqual(record(result, "L", "C")["state"], "always")
        self.assertEqual(record(result, "M", "B")["state"], "partial")
        self.assertEqual(record(result, "M", "D")["state"], "partial")
        self.assertIsNone(record(result, "M", "P1")["component"])

        a = record(result, "E", "A")
        forbidden_a = a["counterfactuals"]["forbid_used"]
        self.assertIsNotNone(forbidden_a)
        self.assertEqual(forbidden_a["trajectory_count"], 0)
        self.assertEqual(forbidden_a["total_cost"], "0")
        self.assertEqual(forbidden_a["witness"], [])
        required_a = a["counterfactuals"]["require_used"]
        self.assertEqual(required_a["trajectory_count"], 1)
        self.assertEqual(required_a["total_cost"], "1/100")

        b = record(result, "M", "B")
        for side in ("require_used", "forbid_used"):
            plan = b["counterfactuals"][side]
            self.assertEqual(plan["trajectory_count"], 1)
            self.assertEqual(plan["total_cost"], "1/100")
            used = identity_keys(plan["witness"][0])
            self.assertEqual(len(used), 3)
            self.assertEqual(len(set(used)), 3)
        self.assertIn(("M", "B"), identity_keys(b["counterfactuals"]["require_used"]["witness"][0]))
        self.assertNotIn(("M", "B"), identity_keys(b["counterfactuals"]["forbid_used"]["witness"][0]))
        self.assertIn(("M", "D"), identity_keys(b["counterfactuals"]["forbid_used"]["witness"][0]))

        static = result["static_photometry"]
        by_exposure = {item["exposure"]: item for item in static["exposures"]}
        self.assertEqual(by_exposure["E"]["detections"], [])
        self.assertEqual(by_exposure["L"]["detections"], [])
        for eid in ("E", "L"):
            cal = by_exposure[eid]["calibration"]
            self.assertEqual(cal["status"], "insufficient_anchors")
            self.assertIsNone(cal["zero_point"])

        middle = by_exposure["M"]
        cal = middle["calibration"]
        self.assertEqual(cal["status"], "insufficient_anchors")
        self.assertIsNone(cal["zero_point"])
        self.assertEqual(
            sorted((a["observation"], a["reference"]) for a in cal["anchors"]),
            [("P1", "S1"), ("P2", "S2")],
        )
        detections = {d["observation"]: d for d in middle["detections"]}
        self.assertEqual(detections["B"]["selected"], "S3")
        self.assertIsNone(detections["B"]["corrected_magnitude"])
        self.assertEqual(detections["D"]["status"], "unmatched")
        self.assertIsNone(detections["F"]["corrected_magnitude"])
        self.assertNotIn("A", by_exposure["E"])

        reordered = copy.deepcopy(data)
        reordered["epochs"].reverse()
        reordered["positions"].reverse()
        reordered["photometry"]["reference"].reverse()
        reordered["photometry"]["exposures"].reverse()
        for exposure in reordered["photometry"]["exposures"]:
            exposure["observations"].reverse()
        self.assertEqual(reduce_moving_cadence(reordered), result)

    def test_middle_epoch_shift_makes_d_unique_and_b_anchor(self):
        data = first_cadence_data("60001.1")
        result = reduce_moving_cadence(data)

        self.assertEqual(record(result, "M", "D")["state"], "always")
        self.assertEqual(record(result, "M", "B")["state"], "never")
        self.assertTrue(record(result, "M", "B")["in_candidates"])
        component = result["components"][0]
        self.assertEqual(component["unconstrained_best"]["trajectory_count"], 1)
        self.assertEqual(component["unconstrained_best"]["total_cost"], "0")

        b = record(result, "M", "B")
        required_b = b["counterfactuals"]["require_used"]
        self.assertEqual(required_b["trajectory_count"], 1)
        self.assertEqual(required_b["total_cost"], "1/25")
        forbidden_b = b["counterfactuals"]["forbid_used"]
        self.assertEqual(forbidden_b["total_cost"], "0")

        by_exposure = {
            item["exposure"]: item
            for item in result["static_photometry"]["exposures"]
        }
        middle = by_exposure["M"]
        self.assertNotIn("D", {d["observation"] for d in middle["detections"]})
        cal = middle["calibration"]
        self.assertEqual(cal["status"], "calibrated")
        self.assertEqual(cal["zero_point"], "2")
        self.assertEqual(
            sorted(a["observation"] for a in cal["anchors"]), ["B", "P1", "P2"]
        )
        detections = {d["observation"]: d for d in middle["detections"]}
        self.assertEqual(detections["B"]["selected"], "S3")
        self.assertEqual(detections["B"]["corrected_magnitude"], "15")
        self.assertEqual(detections["F"]["corrected_magnitude"], "15")

    def test_non_greedy_component_and_static_reassignment(self):
        photometry = {
            "reference": [ref("R", "0", calibrator=False)],
            "threshold": "1",
            "anchor_tolerance": "0.5",
            "exposures": [
                {
                    "id": "E",
                    "observations": [obs("E0", ra="10"), obs("E2", ra="11")],
                },
                {
                    "id": "M",
                    "observations": [
                        obs("M4", ra="12"),
                        obs("M505", ra="13"),
                        obs("M1005", ra="14"),
                    ],
                },
                {
                    "id": "L",
                    "observations": [
                        obs("L8", ra="0.0001"),
                        obs("L20", ra="15"),
                        obs("Q", ra="0.0005"),
                    ],
                },
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "60000"},
                {"exposure": "M", "mjd": "60001"},
                {"exposure": "L", "mjd": "60002"},
            ],
            "positions": [
                position("E", "E0", "0"),
                position("E", "E2", "2"),
                position("M", "M4", "4"),
                position("M", "M505", "5.05"),
                position("M", "M1005", "10.05"),
                position("L", "L8", "8"),
                position("L", "L20", "20"),
                position("L", "Q", "100", "100"),
            ],
            "min_speed": "1",
            "max_speed": "20",
            "middle_tolerance": "0.1",
        }

        before = process_exposures(photometry)
        before_late = next(e for e in before["exposures"] if e["exposure"] == "L")
        before_det = {d["observation"]: d for d in before_late["detections"]}
        self.assertEqual(before_det["L8"]["selected"], "R")
        self.assertIsNone(before_det["Q"]["selected"])

        result = reduce_moving_cadence(data)
        self.assertEqual(len(result["components"]), 1)
        best = result["components"][0]["unconstrained_best"]
        self.assertEqual(best["trajectory_count"], 2)
        self.assertEqual(best["total_cost"], "1/200")
        self.assertEqual(len(best["witness"]), 2)

        always = {"E2", "M505", "M1005", "L8", "L20", "E0"}
        for oid in always:
            self.assertEqual(record(result, oid[0], oid)["state"], "always", oid)
        self.assertEqual(record(result, "M", "M4")["state"], "never")

        m4 = record(result, "M", "M4")
        forced = m4["counterfactuals"]["require_used"]
        self.assertEqual(forced["trajectory_count"], 1)
        self.assertEqual(forced["total_cost"], "0")
        used = [identity_keys(t) for t in forced["witness"]]
        flat = [item for triple in used for item in triple]
        self.assertEqual(len(flat), len(set(flat)))
        self.assertIn(("M", "M4"), flat)

        after_late = next(
            e
            for e in result["static_photometry"]["exposures"]
            if e["exposure"] == "L"
        )
        after_det = {d["observation"]: d for d in after_late["detections"]}
        self.assertNotIn("L8", after_det)
        self.assertNotIn("L20", after_det)
        self.assertEqual(after_det["Q"]["selected"], "R")
        self.assertEqual(after_det["Q"]["status"], "matched")

    def test_sixty_independent_components_are_not_cartesian_expanded(self):
        exposures = [
            {"id": "E", "observations": []},
            {"id": "M", "observations": []},
            {"id": "L", "observations": []},
        ]
        positions = []
        for i in range(60):
            base = 1000 * i
            eid = f"E{i}"
            b = f"B{i}"
            d = f"D{i}"
            cid = f"L{i}"
            exposures[0]["observations"].append(obs(eid))
            exposures[1]["observations"].extend([obs(b), obs(d)])
            exposures[2]["observations"].append(obs(cid))
            positions.extend(
                [
                    position("E", eid, base),
                    position("M", b, f"{base}.9"),
                    position("M", d, f"{base + 1}.1"),
                    position("L", cid, base + 2),
                ]
            )

        data = {
            "photometry": {
                "reference": [ref("S", "180", calibrator=False)],
                "threshold": "0.000001",
                "anchor_tolerance": "0.1",
                "exposures": exposures,
            },
            "epochs": [
                {"exposure": "L", "mjd": "60002"},
                {"exposure": "E", "mjd": "60000"},
                {"exposure": "M", "mjd": "60001"},
            ],
            "positions": list(reversed(positions)),
            "min_speed": "0.75",
            "max_speed": "1.25",
            "middle_tolerance": "0.2",
        }
        result = reduce_moving_cadence(data)
        json.dumps(result)
        self.assertEqual(len(result["components"]), 60)
        self.assertEqual(len(result["candidates"]), 120)
        always = partial = 0
        for component in result["components"]:
            self.assertEqual(component["detection_count"], 4)
            self.assertEqual(component["unconstrained_best"]["trajectory_count"], 1)
            self.assertEqual(component["unconstrained_best"]["total_cost"], "1/100")
            self.assertEqual(len(component["unconstrained_best"]["witness"]), 1)
            for item in component["detections"]:
                state = record(result, item["exposure"], item["observation"])["state"]
                always += state == "always"
                partial += state == "partial"
        self.assertEqual((always, partial), (120, 120))

        shuffled = copy.deepcopy(data)
        shuffled["epochs"] = list(reversed(shuffled["epochs"]))
        shuffled["positions"] = shuffled["positions"][1:] + shuffled["positions"][:1]
        self.assertEqual(reduce_moving_cadence(shuffled), result)

    def test_component_limit_and_identity_validation(self):
        def oversized(middle_count):
            return {
                "photometry": {
                    "reference": [ref("S", "10")],
                    "threshold": "1",
                    "anchor_tolerance": "0.1",
                    "exposures": [
                        {"id": "E", "observations": [obs("A")]},
                        {
                            "id": "M",
                            "observations": [obs(f"M{i}") for i in range(middle_count)],
                        },
                        {"id": "L", "observations": [obs("C")]},
                    ],
                },
                "epochs": [
                    {"exposure": "E", "mjd": "60000"},
                    {"exposure": "M", "mjd": "60001"},
                    {"exposure": "L", "mjd": "60002"},
                ],
                "positions": [
                    position("E", "A", "0"),
                    *[position("M", f"M{i}", "1") for i in range(middle_count)],
                    position("L", "C", "2"),
                ],
                "min_speed": "1",
                "max_speed": "20",
                "middle_tolerance": "0",
            }

        self.assertEqual(len(reduce_moving_cadence(oversized(13))["components"]), 1)
        with self.assertRaises(ValueError):
            reduce_moving_cadence(oversized(14))

        data = first_cadence_data()
        duplicate = copy.deepcopy(data)
        duplicate["positions"][1] = copy.deepcopy(duplicate["positions"][0])
        with self.assertRaises(ValueError):
            reduce_moving_cadence(duplicate)

        missing_position = copy.deepcopy(data)
        missing_position["positions"].pop()
        with self.assertRaises(ValueError):
            reduce_moving_cadence(missing_position)

        unknown_position = copy.deepcopy(data)
        unknown_position["positions"][0]["observation"] = "NOPE"
        with self.assertRaises(ValueError):
            reduce_moving_cadence(unknown_position)

        duplicate_epoch = copy.deepcopy(data)
        duplicate_epoch["epochs"][2] = copy.deepcopy(duplicate_epoch["epochs"][0])
        with self.assertRaises(ValueError):
            reduce_moving_cadence(duplicate_epoch)

        bad_speed = copy.deepcopy(data)
        bad_speed["min_speed"] = "1.3"
        with self.assertRaises(ValueError):
            reduce_moving_cadence(bad_speed)

    def test_strict_fields_ranges_and_photometry_errors(self):
        def assert_invalid(mutate):
            bad = first_cadence_data()
            mutate(bad)
            with self.assertRaises(ValueError):
                reduce_moving_cadence(bad)

        assert_invalid(lambda d: d.update(unexpected=True))
        assert_invalid(lambda d: d.pop("positions"))
        assert_invalid(lambda d: d.update(min_speed="0"))
        assert_invalid(lambda d: d.update(min_speed="1.3"))
        assert_invalid(lambda d: d.update(middle_tolerance="-0.1"))
        assert_invalid(lambda d: d.update(max_speed="1.2345678"))
        assert_invalid(lambda d: d["epochs"][0].update(extra="x"))
        assert_invalid(lambda d: d["epochs"][0].update(mjd="-1"))
        assert_invalid(lambda d: d["epochs"][2].update(mjd="60001"))
        assert_invalid(lambda d: d["positions"][0].update(extra="x"))
        assert_invalid(lambda d: d["positions"][0].update(x="1.2345678"))
        assert_invalid(lambda d: d["photometry"]["reference"][0].pop("magnitude"))
        assert_invalid(lambda d: d["photometry"]["exposures"].pop())
        assert_invalid(lambda d: d.update(positions={}))

        too_many_refs = first_cadence_data()
        too_many_refs["photometry"]["reference"] = [
            ref(f"S{i}", f"{i % 360}") for i in range(257)
        ]
        with self.assertRaises(ValueError):
            reduce_moving_cadence(too_many_refs)

    def test_same_observation_name_across_exposures_and_no_candidate_static_parity(self):
        photometry = {
            "reference": [ref("S", "20")],
            "threshold": "1",
            "anchor_tolerance": "0.2",
            "exposures": [
                {"id": "E", "observations": [obs("same", ra="10")]},
                {"id": "M", "observations": [obs("same", ra="20")]},
                {"id": "L", "observations": [obs("same", ra="30")]},
            ],
        }
        data = {
            "photometry": photometry,
            "epochs": [
                {"exposure": "E", "mjd": "60000"},
                {"exposure": "M", "mjd": "60001"},
                {"exposure": "L", "mjd": "60002"},
            ],
            "positions": [
                position("E", "same", "0", "100"),
                position("M", "same", "0", "-100"),
                position("L", "same", "0", "300"),
            ],
            "min_speed": "0.75",
            "max_speed": "1.25",
            "middle_tolerance": "0.2",
        }
        result = reduce_moving_cadence(data)
        self.assertEqual(len(result["detections"]), 3)
        self.assertTrue(all(item["status"] == "not_in_candidate" for item in result["detections"]))
        self.assertEqual(result["static_photometry"], process_exposures(photometry))


if __name__ == "__main__":
    unittest.main()
