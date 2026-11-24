# data-qd-206 · 天文源目录与不确定性匹配

将观测源按球面角距离和各向同性角误差与参考目录匹配，保留唯一、歧义和无匹配候选；不做图像标定。另提供多次曝光整批接口 `process_exposures`：曝光内一对一关联（最大配对、最小角距代价）、基于定标锚点的逐曝光测光零点与校正星等。

新增离线 `solve_astrometry(data)`：从显式切点和控制星建立 gnomonic 投影上的六参数二维仿射解，执行确定性稳健共识、歧义识别和最终内点凸包责任判定，再转换科学检测 RA/Dec。详见 [二维仿射天体测量契约](docs/astrometry.md)。

另提供仅依赖标准库的离线跟拍工作簿 `FollowUpWorkbook` / `create_workbook` / `open_workbook`：持久化测光批次、整批质量更正、精确符号检验与 BH 校正、按观测组优先序维护带永久令牌的待执行指令，并登记不可撤销、可幂等重试的实拍记录。详见 [离线跟拍工作簿契约](docs/followup.md)。

新增离线测光发布件 `build_photometric_release(data)` / `assess_fade(release, query)`：以几何最优方案签发可信实测、核验清洁的未检出灵敏度证明、按转移连通分量枚举零点故障假设（不删或恰好删一条），并仅用发布件在新进程里给出精确有理、可复核、可追溯检测 ID 的变暗下界与确认结论。发布件为版本化 JSON，不含原始曝光行、坐标或候选边表。详见 [测光发布件契约](docs/release.md)。

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
.venv/bin/python -I examples/run_followup.py
.venv/bin/python -I examples/run_release.py
```

上述命令已在本地缓存和独立环境真实执行；样例输出写入 `run/sample.json`、`run/exposures.json` 与 `run/photometric_release.json`，已有目标拒绝覆盖。若暂不安装，也可在仓库根目录运行 `python3 -I run_tests.py` 执行全部标准库测试，或直接运行 `python3 -I examples/run_followup.py`。天体测量离线示例可运行 `python3 -I examples/run_astrometry.py`（脚本自行注入 `src/`，输出 `run/astrometry.json`）。API可这样调用：

```python
from star_crossmatch import process, process_exposures, solve_astrometry
from star_crossmatch import build_photometric_release, assess_fade
from star_crossmatch.io import load, publish
publish("run/another.json", process(load("fixtures/input.json")))

# 发布件可 json.dumps 后交接;交接方在新进程仅凭 JSON 即可 assess_fade
release = build_photometric_release(batch_with_epochs_transfers_limits)
verdict = assess_fade(json.loads(json.dumps(release)), fade_query)

# astrometry_input 的完整字段见 docs/astrometry.md
astrometry_result = solve_astrometry(astrometry_input)
```

跟拍工作簿使用本地 JSON 持久化，关闭后可用 `open_workbook(path)` 重开：

```python
from star_crossmatch import create_workbook, open_workbook

wb = create_workbook("run/followup.json", batch)
print(wb.candidate_results())       # 各源 n/k/p/q 与候选/证据不足状态
print(wb.pending_instructions())    # 带永久令牌和发令证据的待执行指令
receipt = wb.record_execution(token, "shot-record-id")
wb.update_quality_flags([{"measurement_id": "m1", "usable": False}])
wb.close()
wb = open_workbook("run/followup.json")
```

完整字段、精确有理数口径、撤销/新发规则和冲突身份处理见
[docs/followup.md](docs/followup.md)。

样例摘要：

```json
{"ambiguous": 0, "matched": 1, "observations": 2, "reference": 2, "unmatched": 1}
```

现有匹配、多次曝光与离线跟拍业务测试均通过（含完整容量 2 验收剧本、精确有理数边界、重开持久化和冲突身份），并核对非法输入拒绝；样例数据全部原创合成，见 [来源](fixtures/SOURCES.md)。

## 边界

业务结果仅覆盖契约中定义的有限模型。后续需求可沿已有身份和证据字段扩展，不把当前样例当成行业标准或生产认证。
