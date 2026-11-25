#!/usr/bin/env python3
"""离线距离后验示例：含负/正视差、共享零点与源对排序，不访问网络。

在仓库根目录运行::

    python3 -I examples/run_distance_catalog.py

结果通过 io.publish 原子新建到 run/distance_catalog.json；若目标已存在会拒绝覆盖。
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from star_crossmatch import infer_distance_catalog  # noqa: E402
from star_crossmatch.io import publish  # noqa: E402

DATA = {
    "scale_length_pc": "1200",
    "zero_point": {"mean_mas": "0.08", "sigma_mas": "0.15"},
    "sources": [
        # 负视差源：不截断、不取倒数，后验由先验与似然共同决定，距离严格为正。
        {"id": "A", "parallax_mas": "-0.20", "sigma_mas": "0.50"},
        {"id": "B", "parallax_mas": "0.75", "sigma_mas": "0.35"},
        {"id": "C", "parallax_mas": "0.10", "sigma_mas": "0.42"},
    ],
    "pairs": [
        {"id": "A-before-B", "nearer": "A", "farther": "B"},
        {"id": "C-before-B", "nearer": "C", "farther": "B"},
        {"id": "B-before-A", "nearer": "B", "farther": "A"},
    ],
}


def main():
    result = infer_distance_catalog(DATA)
    print("format:", result["format"])
    print("zero point posterior:", result["zero_point"])
    for row in result["sources"]:
        print(
            row["id"],
            f"median={row['median_pc']:.3f} pc",
            f"68.27% CI=({row['ci_lower_pc']:.3f},"
            f"{row['ci_upper_pc']:.3f}) pc",
        )
    for pair in result["pairs"]:
        print(pair["id"], "P(r_near<r_far) =", pair["probability"])
    target = ROOT / "run" / "distance_catalog.json"
    publish(target, result)
    print(f"written {target}")


if __name__ == "__main__":
    main()
