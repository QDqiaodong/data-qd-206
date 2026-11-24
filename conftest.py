"""Pytest 本地开发配置：使用 src/ 布局时免安装导入包。

正式运行只依赖标准库；不使用 pytest 时可执行 ``python3 -I run_tests.py``。
"""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
