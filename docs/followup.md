# 离线跟拍工作簿契约

跟拍工作簿服务于一个有限跟拍窗口：从一批测光数据计算亮源符号检验和
Benjamini–Hochberg（BH）校正，按观测组预先给定的优先序维护待执行指令，
并登记已经实拍、不可撤销的指令。它不做坐标关联、零点估计或望远镜通信。

工作簿只依赖 Python 标准库。所有成功修改都会原子写入同一个本地 JSON 文件；
不要求、也不应让多个进程同时写同一个工作簿。

## 初始化与重开

```python
from star_crossmatch import FollowUpWorkbook, create_workbook, open_workbook

batch = {
    "sources": [
        {"id": "X", "baseline_magnitude": "15.0", "priority": 1},
        # ...
    ],
    "measurements": [
        {
            "id": "x0",
            "source_id": "X",
            "mjd": "59000.0",
            "magnitude": "14.9",
            "sigma": "0.01",
            "usable": True,
        },
        # ...
    ],
    "min_valid_times": 5,
    "fdr_threshold": "0.04",
    "capacity": 2,
}

# 两种写法等价：
wb = create_workbook("run/window.json", batch)
wb = FollowUpWorkbook.create("run/window.json", batch)

# 也可以用具名参数：
wb = FollowUpWorkbook.create(
    "run/window.json",
    sources=batch["sources"],
    measurements=batch["measurements"],
    min_valid_times=5,
    fdr_threshold="0.04",
    capacity=2,
)

wb.close()                  # 所有成功修改此前均已落盘
wb = open_workbook("run/window.json")
```

路径直接表示工作簿 JSON 文件，是否使用 `.json` 后缀均可。若路径是已存在目录，
工作簿文件为该目录下的 `followup_workbook.json`。初始化拒绝覆盖已存在的工作簿文件。

ID 使用仓库既有语法 `[A-Za-z0-9_:-]{1,64}`；十进制输入使用既有语法
（最多 9 位整数、6 位小数）。初始化边界：

- `sources`：1..1000 项；`id` 唯一；`priority` 必须是 1..源数的全批唯一排列。
- `measurements`：1..10000 项；`id` 全批唯一；`source_id` 必须已登记。
- `mjd` 非负；`sigma` 为正；`usable` 必须是 JSON 布尔值。
- 每个源最多有 128 个按精确有理数判定的不同输入时刻。
- `min_valid_times`：整数 3..128。
- `fdr_threshold`：十进制字符串，精确范围 `(0, 1]`。
- `capacity`：整数 1..源数。

无效初始化抛 `ValueError`（具体类型可能是 `FollowUpError`），且不会创建工作簿
状态文件。字段别名（如 `baseline`、`source`、`record_id`/`shot_id`）只是为了方便
调用；上文名称是规范名称。

## 精确测光统计

同一源、同一数值 MJD 的可用记录按逆方差权重合并：

\[
\bar m = \frac{\sum_i m_i/\sigma_i^2}{\sum_i 1/\sigma_i^2}.
\]

MJD 相等、星等/标准差解析、逆方差合并以及与基准星等的比较全部使用
`fractions.Fraction`。因此 `"1.0"` 与 `"1.000000"` 是同一时刻；合并星等严格小于
基准记为偏亮，严格大于记为未偏亮，恰好相等时该时刻不计入 `n` 或 `k`。

非相等时刻数为 `n`，偏亮时刻数为 `k`。当 `n < min_valid_times` 时状态为
`insufficient_evidence`，`p`、`q` 均为 `null`，不参加全批 BH 校正。其余源的单侧
符号检验 p 值为公平硬币独立抛 `n` 次、至少 `k` 次正面的精确概率：

\[
p = 2^{-n}\sum_{j=k}^{n}\binom{n}{j}.
\]

BH 校正在当前全部合格源上进行。先按 `(p, source_id)` 排序并从 1 开始排名，再取
从该排名到末尾的 `合格源数 * p_j / j` 的最小值，上限为 1。p 值相等的源虽然按 ID
决定排名，但反向累计最小值保证它们 q 值相同。只有精确满足 `q <= fdr_threshold`
的源状态为 `candidate`。

### 读取当前结果

```python
rows = wb.candidate_results()   # 全部源，按 priority 排序
rows["B"]                       # IndexedList 也可用源 ID 查找
wb.current_candidates()         # 仅 q <= threshold 的源
wb.get_source("B")
```

每行形如：

```python
{
    "source_id": "B",
    "source": "B",
    "id": "B",                 # source_id 的别名
    "priority": 2,
    "n": 5,
    "k": 5,
    "p": "1/32",               # 约分分数；整数分母时输出整数字符串
    "q": "5/128",
    "candidate": True,
    "is_candidate": True,
    "insufficient_evidence": False,
    "status": "candidate",      # candidate | not_candidate | insufficient_evidence
}
```

当前统计每次从测光标记重新计算；它与指令中保存的“发令时证据”是两套查询结果。

## 待执行指令与整批更正

```python
pending = wb.pending_instructions()
# 也可写 wb.instructions() / wb.pending() / wb.read()["instructions"]
```

指令形如：

```python
{
    "token": "fu-000000000002",  # 单调、持久、永不复用
    "source_id": "B",
    "source": "B",                # source_id 的别名
    "n": 5,
    "k": 5,
    "p": "1/32",
    "q": "5/128",
}
```

窗口剩余名额为 `capacity - len(executed_records())`。工作簿从尚未执行的当前候选中，
按观测组 `priority` 升序取前“剩余名额”个源作为计划。已在计划内且仍属于该集合的源
保留原令牌和原发令证据；离开集合的指令立即失效；以后重新进入时一定获得新令牌。

整批质量更正只按测量 ID 设置 `usable`：

```python
result = wb.update_quality_flags([
    {"measurement_id": "a0", "usable": False},
    {"measurement_id": "a1", "usable": False},
])

result["revoked"]  # 本次失效的原指令
result["issued"]   # 本次新发的指令
result["new"]      # result["issued"] 的别名
```

也接受 `{"a0": False, "a1": False}` 映射。同一批中同一测量 ID 不得出现两次；
只要有未知测量 ID、重复指针或无效布尔值，整批更正抛 `ValueError`/
`FollowUpError`，测光标记和跟拍计划均不变。

## 登记实际拍摄

```python
receipt = wb.record_execution(token, "shot-X-001")
# 兼容关键字：record_execution(token, record_id="...")
# 或 record_execution(token, shot_id="...")
```

成功后：

- 令牌立即从 `pending_instructions()` 消失；
- 原指令和实拍记录 ID 作为一条独立记录保存；
- `execution_count` 增加 1，剩余容量减少 1；
- 同一 `(token, record_id)` 原样重试返回相等的原执行记录，不再次扣容量；
- 同一令牌换另一实拍 ID、同一实拍 ID 配另一令牌、已撤销/未知令牌均抛
  `FollowUpConflict`（`ValueError` 子类）；
- 已执行记录不能删除或改写成未执行。

执行记录形如：

```python
{
    "token": "fu-000000000001",
    "record_id": "shot-X-001",
    "shot_id": "shot-X-001",
    "capture_id": "shot-X-001",
    "instruction": { ... 原发令证据 ... },
}
```

使用 `wb.executed_records()`（别名 `executions()`/`records()`）读取全部执行记录。
`wb.read()` 一次返回候选、当前候选、待执行指令、执行记录和容量状态。

## 原子性与重开

每次成功的更正或登记都在同一个新状态对象上完成全部校验和计划调和，随后写入同目录
临时文件、`fsync`、`os.replace` 并 `fsync` 目录。因此：

- 无效更正/登记不会修改当前内存状态或现有 JSON 文件；
- 不会出现“执行记录已保存但同一指令仍待执行”的中间状态；
- 不会出现“质量标记已改但旧计划未调和”的中间状态；
- 关闭、进程结束或重新调用 `open_workbook(path)` 后，所有成功修改均可见；
- 不提供多进程并发写保证。

## 示例与测试

```sh
python3 -I examples/run_followup.py
python3 -I run_tests.py                 # 免安装，自动加入 src/
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

安装为包后可直接按 README 中的虚拟环境方式运行 `python -I -m unittest ...`。
`tests/test_followup.py` 覆盖题目要求的完整容量 2 验收剧本、重开、幂等回执、
冲突身份、整批失败回滚、精确有理数时刻/合并、行重排、p 并列、q 恰等于阈值以及
发令证据保留。
