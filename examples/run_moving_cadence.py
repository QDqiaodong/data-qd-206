#!/usr/bin/env python3
"""离线最小示例:三帧移动轨迹审计与静态测光交接。

在仓库根目录运行::

    python3 -I examples/run_moving_cadence.py

结果通过 ``io.publish`` 原子新建到 ``run/moving_cadence.json``;目标已存在
时会拒绝覆盖。示例全部使用合成坐标,不访问网络。
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from star_crossmatch import reduce_moving_cadence  # noqa: E402
from star_crossmatch.io import publish  # noqa: E402


def observation(oid, ra):
    return {
        "id": oid,
        "ra": ra,
        "dec": "0",
        "sigma": "0.001",
        "instrumental_magnitude": "13",
    }


DATA = {
    "photometry": {
        "reference": [
            {"id": "S1", "ra": "10", "dec": "0", "magnitude": "15", "calibrator": True},
            {"id": "S2", "ra": "20", "dec": "0", "magnitude": "15", "calibrator": True},
            {"id": "S3", "ra": "30", "dec": "0", "magnitude": "15", "calibrator": True},
        ],
        "threshold": "1",
        "anchor_tolerance": "0.2",
        "exposures": [
            {"id": "E", "observations": [observation("A", "0")]},
            {
                "id": "M",
                "observations": [
                    observation("B", "30"),
                    observation("D", "40"),
                    observation("P1", "10"),
                    observation("P2", "20"),
                ],
            },
            {"id": "L", "observations": [observation("C", "1")]},
        ],
    },
    "epochs": [
        {"exposure": "E", "mjd": "60000"},
        {"exposure": "M", "mjd": "60001"},
        {"exposure": "L", "mjd": "60002"},
    ],
    "positions": [
        {"exposure": "E", "observation": "A", "x": "0", "y": "0"},
        {"exposure": "M", "observation": "B", "x": "0.9", "y": "0"},
        {"exposure": "M", "observation": "D", "x": "1.1", "y": "0"},
        {"exposure": "M", "observation": "P1", "x": "10", "y": "10"},
        {"exposure": "M", "observation": "P2", "x": "11", "y": "10"},
        {"exposure": "L", "observation": "C", "x": "2", "y": "0"},
    ],
    "min_speed": "0.75",
    "max_speed": "1.25",
    "middle_tolerance": "0.2",
}


if __name__ == "__main__":
    result = reduce_moving_cadence(DATA)
    best = result["global_best"]
    print("global best:", best["trajectory_count"], best["total_cost"])
    for detection in result["detections"]:
        print(
            detection["exposure"],
            detection["observation"],
            detection["state"],
            detection["component"],
        )
    middle = next(
        exposure
        for exposure in result["static_photometry"]["exposures"]
        if exposure["exposure"] == "M"
    )
    print("middle calibration:", middle["calibration"]["status"])
    publish("run/moving_cadence.json", result)
    print("written run/moving_cadence.json")
