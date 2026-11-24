#!/usr/bin/env python3
"""离线 astrometry 示例：合成切点、控制星和科学检测，不访问网络或图像。

在仓库根目录运行::

    python3 -I examples/run_astrometry.py

结果通过 io.publish 原子新建到 run/astrometry.json；若目标已存在会拒绝覆盖。
"""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from star_crossmatch import solve_astrometry  # noqa: E402
from star_crossmatch.io import publish  # noqa: E402

RA0 = 0.0
DEC0 = 0.0
# xi  = A*x + B*y + C；eta = D*x + E*y + F，单位为度。
AFFINE = (0.010, 0.002, 0.050, 0.001, 0.012, -0.030)


def inverse_gnomonic(xi_deg, eta_deg):
    xi = math.radians(xi_deg)
    eta = math.radians(eta_deg)
    dec0 = math.radians(DEC0)
    rho = math.hypot(xi, eta)
    if rho == 0.0:
        return RA0, DEC0
    c = math.atan(rho)
    dec = math.asin(
        math.cos(c) * math.sin(dec0)
        + eta * math.sin(c) * math.cos(dec0) / rho
    )
    ra = math.radians(RA0) + math.atan2(
        xi * math.sin(c),
        rho * math.cos(dec0) * math.cos(c)
        - eta * math.sin(dec0) * math.sin(c),
    )
    return math.degrees(ra) % 360.0, math.degrees(dec)


def sky(x, y):
    a, b, c, d, e, f = AFFINE
    return inverse_gnomonic(a * x + b * y + c, d * x + e * y + f)


def control(cid, x, y, *, bad=False):
    if bad:
        ra, dec = 12.345, -6.789
    else:
        ra, dec = sky(x, y)
    return {
        "id": cid,
        "x": f"{x:.6f}",
        "y": f"{y:.6f}",
        "ra": f"{ra:.6f}",
        "dec": f"{dec:.6f}",
    }


DATA = {
    "tangent_point": {"ra": "0", "dec": "0"},
    "controls": [
        control("c00", 0, 0),
        control("c10", 10, 0),
        control("c20", 20, 0),
        control("c01", 0, 10),
        control("c11", 10, 10),
        control("c21", 20, 10),
        control("c02", 0, 20),
        control("c12", 10, 20),
        control("c22", 20, 20),
        # 错配/坏测量：共识筛选应将其标为非内点，且不让它扩大凸包。
        control("bad-1", 25, 25, bad=True),
    ],
    "detections": [
        {"id": "det-inside", "x": "10.000000", "y": "10.000000"},
        {"id": "det-edge", "x": "20.000000", "y": "10.000000"},
        {"id": "det-outside", "x": "40.000000", "y": "40.000000"},
    ],
    "residual_limit_arcsec": "0.05",
    "min_inliers": 6,
}


def main():
    result = solve_astrometry(DATA)
    print("status:", result["status"])
    print("inliers:", result["inlier_count"])
    for det in result["detections"]:
        print(det["id"], det["status"], det["ra"], det["dec"])
    target = ROOT / "run" / "astrometry.json"
    publish(target, result)
    print(f"written {target}")


if __name__ == "__main__":
    main()
