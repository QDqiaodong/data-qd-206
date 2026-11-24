#!/usr/bin/env python3
"""无需安装即可运行全部测试的标准库入口。

用法::

    python3 -I run_tests.py
    python3 -I run_tests.py -v

已按 README 安装包后，也可以继续使用
``python -I -m unittest discover -s tests -v``。
"""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def main():
    verbosity = 2 if "-v" in sys.argv[1:] or "--verbose" in sys.argv[1:] else 1
    suite = unittest.defaultTestLoader.discover(
        str(ROOT / "tests"), pattern="test_*.py"
    )
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
