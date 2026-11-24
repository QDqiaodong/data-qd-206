"""离线最小演示:测光发布件构建、JSON 交接与变暗复核。

运行::

    python3 -I examples/run_release.py

脚本只依赖标准库,把发布件 JSON 写入 run/photometric_release.json,
并用纯 dict(模拟交接方在新进程里只拿到发布件)完成两次变暗查询。
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from star_crossmatch import assess_fade, build_photometric_release
from star_crossmatch.io import publish


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


def build(cd_delta):
    photometry = {
        "reference": [
            ref("C1", "10", "0", "15", True),
            ref("C2", "20", "0", "15", True),
            ref("C3", "30", "0", "15", True),
            ref("S", "50", "0", "18"),
        ],
        "threshold": "1",
        "anchor_tolerance": "0.5",
        "exposures": [
            {"id": "A", "observations": [
                obs("k1", "10", "0", "0.001", "13"),
                obs("k2", "20", "0", "0.001", "13"),
                obs("k3", "30", "0", "0.001", "13"),
            ]},
            {"id": "B", "observations": [obs("b1", "80", "0", "0.001", "13")]},
            {"id": "C", "observations": [obs("sC", "50", "0", "0.001", "16")]},
            {"id": "D", "observations": [obs("d1", "70", "0", "0.001", "13")]},
        ],
    }
    return {
        "photometry": photometry,
        "epochs": [
            {"exposure": "A", "mjd": "1"},
            {"exposure": "B", "mjd": "2"},
            {"exposure": "C", "mjd": "3"},
            {"exposure": "D", "mjd": "4"},
        ],
        "transfers": [
            transfer("t1", "A", "B", "1"),
            transfer("t2", "B", "C", "1"),
            transfer("t3", "A", "C", "3"),
            transfer("t4", "C", "D", cd_delta),
        ],
        "limits": [
            {"exposure": "D", "reference": "S",
             "limiting_instrumental_magnitude": "16"}
        ],
    }


def main():
    query = {
        "reference": "S",
        "bright_exposure": "C",
        "later_exposure": "D",
        "minimum_fade": "1",
    }
    for label, cd_delta in (("自洽主剧本(C->D=1)", "1"), ("冲突变体(C->D=0.5)", "0.5")):
        release = build_photometric_release(build(cd_delta))
        print(f"== {label} ==")
        for row in release["exposures"]:
            print(
                " ", row["exposure"], row["zero_point_status"],
                "zero_point=", row["zero_point"],
            )
        print(
            "  可行假设:",
            [h["excluded_transfer"] for h in release["components"][0]["hypotheses"]],
        )
        # 交接方视角:只拿经 JSON 往返的纯 dict,没有任何原始输入或进程缓存。
        handed = json.loads(json.dumps(release))
        verdict = assess_fade(handed, query)
        print("  变暗结论:", verdict["status"], verdict)
        if label.startswith("自洽"):
            publish("run/photometric_release.json", release)
    print("written run/photometric_release.json")


if __name__ == "__main__":
    main()
