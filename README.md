# data-qd-206 · 天文源目录与不确定性匹配

将观测源按球面角距离和各向同性角误差与参考目录匹配，保留唯一、歧义和无匹配候选；不做图像标定。另提供多次曝光整批接口 `process_exposures`：曝光内一对一关联（最大配对、最小角距代价）、基于定标锚点的逐曝光测光零点与校正星等。

本仓主要交付是Python业务API和可供下游使用的结构化数据，不是通用CLI、报表或CRUD模板。字段、身份、单位、上限与失败边界见 [数据契约](docs/contract.md)。

## 安装与验证

Python >=3.11，实测CPython3.12.14、uv0.12.19。运行依赖仅标准库，构建固定setuptools80.9.0、wheel0.45.1；每仓独立 `.venv`。

```sh
uv venv --offline --python 3.12 .venv
uv pip install --offline --python .venv/bin/python -r requirements-build.lock
uv pip install --offline --python .venv/bin/python --no-build-isolation --no-deps .
uv pip check --python .venv/bin/python
.venv/bin/python -I -m unittest discover -s tests -v
.venv/bin/python -I examples/run_sample.py
.venv/bin/python -I examples/run_exposures.py
```

上述命令已在本地缓存和独立环境真实执行；样例输出写入 `run/sample.json` 与 `run/exposures.json`，已有目标拒绝覆盖。API可这样调用：

```python
from star_crossmatch import process, process_exposures
from star_crossmatch.io import load, publish
publish("run/another.json", process(load("fixtures/input.json")))
```

样例摘要：

```json
{"ambiguous": 0, "matched": 1, "observations": 2, "reference": 2, "unmatched": 1}
```

已通过 20 项业务测试（含多次曝光关联、定标零点、并列方案证据与暴力枚举对照），并核对非法输入拒绝；样例数据全部原创合成，见 [来源](fixtures/SOURCES.md)。

## 边界

业务结果仅覆盖契约中定义的有限模型。后续需求可沿已有身份和证据字段扩展，不把当前样例当成行业标准或生产认证。
