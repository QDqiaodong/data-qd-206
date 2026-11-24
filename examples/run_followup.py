#!/usr/bin/env python3
"""离线跟拍工作簿示例：精确候选、发令/撤销、幂等登记与重开。

直接运行：

    python -I examples/run_followup.py

脚本只使用标准库；它会在系统临时目录中创建一个真实 JSON 工作簿，并在结束时
打印路径。每次运行都使用新的临时文件，不会覆盖仓库中的样例。
"""
from pathlib import Path
import sys
import tempfile

# 允许不安装包时直接从仓库源码运行。
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from star_crossmatch.followup import create_workbook, open_workbook  # noqa: E402


def batch():
    priorities = {"X": 1, "B": 2, "D": 3, "A": 4, "E": 5}
    counts = {
        "X": (8, "bright"),
        "B": (5, "bright"),
        "D": (7, "bright"),
        "A": (6, "bright"),
        "E": (5, "faint"),
    }
    sources = [
        {"id": sid, "baseline_magnitude": "15", "priority": priority}
        for sid, priority in priorities.items()
    ]
    measurements = []
    for sid, (count, kind) in counts.items():
        for index in range(count):
            measurements.append(
                {
                    "id": f"{sid.lower()}{index}",
                    "source_id": sid,
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


def main():
    directory = tempfile.mkdtemp(prefix="star-followup-example-")
    path = Path(directory) / "followup_workbook.json"

    wb = create_workbook(path, batch())
    initial = wb.candidate_results()
    print("初始 B: n={n} k={k} p={p} q={q} candidate={candidate}".format(**initial["B"]))
    print("初始待执行:", [(i["source_id"], i["token"]) for i in wb.pending_instructions()])

    x_token = wb.pending_instructions()["X"]["token"]
    submission = (x_token, "shot-X-001")
    receipt = wb.record_execution(*submission)
    print("已登记 X:", receipt["record_id"], receipt["token"])
    wb.close()

    # 模拟未收到回执：重开后用同一令牌和同一实拍 ID 原样重试。
    wb = open_workbook(path)
    retried = wb.record_execution(*submission)
    assert retried == receipt
    assert wb.execution_count == 1

    changed = wb.update_quality_flags(
        [
            {"measurement_id": "a0", "usable": False},
            {"measurement_id": "a1", "usable": False},
        ]
    )
    b_row = wb.candidate_results()["B"]
    print("屏蔽 A 两个时刻后 B q =", b_row["q"])
    print(
        "撤销:", [i["source_id"] for i in changed["revoked"]],
        "新发:", [i["source_id"] for i in changed["issued"]],
        "当前待执行:", [i["source_id"] for i in wb.pending_instructions()],
    )
    assert b_row["q"] == "1/24"
    assert [i["source_id"] for i in wb.pending_instructions()] == ["D"]

    wb.close()
    reopened = open_workbook(path)
    assert reopened.record_execution(*submission) == receipt
    print("重开和幂等重试均通过；工作簿保留在", path)


if __name__ == "__main__":
    main()
