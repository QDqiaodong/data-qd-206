# 天文目录匹配契约

使用者是天文观测目录整理人员，输入reference（1..1000项id/ra/dec/magnitude）与observations（1..2000项id/ra/dec/sigma），threshold为正数。各数组id唯一ASCII；数值是最多9位整数和6位小数的十进制字符串。ra取[0,360)度，dec取[-90,90]度，sigma为正的各向同性角误差（度）。magnitude只保留为来源元数据，不参与筛选。

角距离按球面haversine计算，采用Python双精度三角函数；可跨赤经0/360与极区。距离<=threshold*sigma才是候选，边界按实际浮点结果比较，无额外epsilon。distance为浮点值字符串，不声称精确有理角距离。候选按距离、reference id排序；恰一候选才selected，多个候选ambiguous且selected=null，零候选unmatched。各观测独立，一个参考源可被多个观测命中，不进行全局一对一匹配。

这不是误差椭圆、相机重建或天文图像叠加；没有自行运动/历元换算、球面协方差、亮度联合推断。后续历元传播和源关联可接入reference/observation身份与候选证据，不用改造为通用报表。样例仅为合成角坐标。process输入错误抛ValueError；io.load严格拒绝重复JSON键，io.publish原子新建目标，已有文件或链接拒绝覆盖。

# 多次曝光整合契约（process_exposures）

输入对象字段恰为reference、threshold、anchor_tolerance、exposures。reference为1..1000行，每行恰含id/ra/dec/magnitude/calibrator，calibrator为JSON布尔（是否可用于定标）；threshold含义同process（正十进制）；anchor_tolerance为非负十进制（星等）。exposures为1..2000项，每项恰含id与非空observations（1..2000条）；观测行恰含id/ra/dec/sigma/instrumental_magnitude（仪器星等）。参考星ID与曝光ID各自全批唯一，观测ID仅在本曝光内唯一，全批观测合计1..2000条。ID与数值语法、坐标范围、sigma/threshold正值约束与process相同。任何无效字段、类型、身份、坐标或边界值抛ValueError；不修改传入数据，不返回部分成功结果。

每次曝光独立关联，其他曝光可再次使用同一参考星。候选判定与process相同：haversine浮点距离<=threshold*sigma，亮度不参与入围。每条候选边的整数代价=floor(浮点角距离（度）*3600000+0.5)（毫角秒），先按原半径入围再取整。曝光内每条检测与每颗参考星至多参与一对，检测可留空。几何最优：先最大化配对数，再最小化整数代价之和；配对数与总代价都相同的方案为同等最优，不按ID顺序取舍。全部最优方案由最小费用流的零费用残差环流刻画（强连通分量），不枚举方案。

测光零点只取自所有几何最优方案中都存在且参考星calibrator=true的配对（锚点）。每个锚点偏移=参考星等−仪器星等。锚点≥3才计算偏移中位数（偶数取中间两值的均值）；|偏移−中位数|≤anchor_tolerance为支持；支持者≥3且为全体锚点的严格多数时采用该中位数为零点（status=calibrated）。否则区分锚点不足（insufficient_anchors，锚点<3）与一致性不足（no_consensus），零点为null。偏移、容差判断、中位数与星等残差均按十进制输入的精确值计算，不用浮点近似破平局。零点有效时，所有检测（含无候选、未获分配、歧义检测）的校正星等=仪器星等+零点。

零点有效时，仅在几何最优方案中再比较所有已配对检测的|校正星等−参考星等|之和，保留和最小的全部方案；零点无效时保留全部几何最优方案。测光改变的关联不回流参与本次定标。

输出为可JSON序列化对象{"exposures":[...]}，曝光按曝光ID字典序。每曝光：exposure（ID）；calibration={status（calibrated/insufficient_anchors/no_consensus）、zero_point（十进制字符串或null）、anchors}，anchors为用于判断的全部锚点（含未支持中位数的锚点），每项{observation、reference、offset（十进制字符串）}，按（检测ID，参考ID）排序；detections按检测ID字典序。每检测：observation（ID）；candidates为全部几何候选，每项{reference、distance（浮点字符串，同process）、cost（整数毫角秒）}，按（角距，参考ID）排序；possible为最终最优方案中可能配到的参考星ID列表，按ID排序；may_be_unassigned表示是否存在让其未分配的最优方案；selected仅当possible唯一且may_be_unassigned=false时为该参考星ID，否则为null；status为matched（已选中）、ambiguous（最优方案之间有歧义）、unassigned（有几何候选但任何最优方案都不配对）、unmatched（无几何候选）之一；corrected_magnitude为十进制字符串或null。按证据身份比较时，重排输入的参考行、曝光或检测行不改变输出。

# 三帧移动轨迹归约与静态测光交接（reduce_moving_cadence）

新增 `reduce_moving_cadence(data)`：输入恰含 `photometry`、`epochs`、`positions`、`min_speed`、`max_speed`、`middle_tolerance`。`photometry` 复用 `process_exposures` 契约，但额外限定恰好 3 帧、每帧 1..120 检测、1..256 参考源；`epochs` 按实际 MJD 确定早/中/晚，`positions` 给出统一平面角秒坐标。候选由早/中/晚各一条检测组成，精确校验端点速度闭区间和实际 MJD 线性插值容差，代价为中间残差平方的约分有理数字符串。共享检测的候选构成冲突分量，每分量参与候选的不同检测不得超过 15 条，超限整批 `ValueError`。

每个分量先最大化轨迹数、再最小化总代价。检测审计区分每个最佳方案必用、部分方案使用和无最佳方案使用，并给出参与候选检测在“硬要使用/硬要不用”约束下的分量最佳值及完整见证。必用检测从静态测光输入剔除；部分方案使用的检测保留并参与配对/校正，但其配对不得作为零点锚点；其余检测不受限。静态结果在剔除后的输入上重跑原球面候选、整数角距代价、一对一几何最优、锚点/零点和星等破平局规则。完整 JSON 字段、见证格式、排序和示例见 [三帧移动巡天归约契约](moving_cadence.md)。

# 二维仿射天体测量契约（solve_astrometry）

新增 `solve_astrometry(data)`：输入 `tangent_point={ra,dec}`、`controls=[{id,x,y,ra,dec}]`、`detections=[{id,x,y}]`、正十进制字符串 `residual_limit_arcsec` 和整数 `min_inliers`。ID/十进制语法/字段严格性/数组与坐标边界沿用本文件前述约定；控制星与检测 ID 分别唯一。切点必须显式给出，本 API 不从探测器坐标反求指向。

控制星 RA/Dec 先以切点为中心做 gnomonic 投影得到度单位 `xi,eta`，再拟合 `xi=a*x+b*y+c`、`eta=d*x+e*y+f`。确定性共识按 ID 排序，小控制星集穷举非共线三元组，大集合使用文档固定种子的可复现抽样；以最终仿射解的球面角残差 `<=residual_limit_arcsec` 认定内点，全部内点 LS 重拟合至集合稳定。合格解先比内点数、再比内点残差平方和。同分最佳解若对任一检测的球面位置差严格超过残差上限，返回 `status="ambiguous_solution"`；秩不足为 `rank_deficient`，合格内点不足为 `insufficient_inliers`。

有效解仅负责最终内点探测器坐标的精确 Fraction 凸包（含边界）：内部/边界检测返回归一到 `[0,360)` 的 RA/Dec，凸包外检测为 `outside_footprint` 且坐标 null；无有效解时所有检测坐标 null。控制星与检测均按 ID 排序，输出投影坐标、角秒残差、内点状态、仿射参数和检测状态。完整字段、序列化、候选/平局排序和示例见 [二维仿射天体测量契约](astrometry.md)。

# 测光发布件与变暗判据契约（build_photometric_release / assess_fade）

新增 `build_photometric_release(data)` 与 `assess_fade(release, query)`：在 `process_exposures` 的几何最优方案与锚点规则之上，叠加历元、有向转移记录与灵敏度证明，产出版本化、可 JSON 往返、可跨进程离线查询的测光发布件。检测身份只认所有几何最优方案中都锁定的同一条检测（不做星等破平局）；灵敏度证明一旦被任一原始几何候选命中即不可用作未检出证据。转移按无向连通分量建模，自洽分量只保留“不排除”假设，矛盾分量枚举全部可行的单条排除假设，需排除两条及以上者标为无解且不阻塞别组；帧零点仅在分量内每个可行假设给出同一确定绝对值时发布。`assess_fade` 按“观测证据门→定标门→联合假设→精确有理变暗下界”的固定顺序给出可区分、可追溯的结论。输入规模（1..32 曝光、1..256 参考源、全批 1..2000 检测、0..128 转移、0..5000 极限行）、发布件 JSON 结构、状态枚举、并列规则与交接方式见 [测光发布件契约](release.md)。

# 整批视差贝叶斯距离推断契约（infer_distance_catalog）

新增 `infer_distance_catalog(data)`：输入恰含 `sources`（1..64 行 `id/parallax_mas/sigma_mas`，视差可带符号、标准差严格为正）、`zero_point`（`mean_mas/sigma_mas`，后者非负且允许 0）、正十进制 `scale_length_pc`、`pairs`（0..256 行 `id/nearer/farther`，端点为不同的已声明源）。模型为共享加性零点 `z ~ Normal(mean,sigma^2)`、距离先验 `p(r|L)=r^2 exp(-r/L)/(2L^3)`、视差 `y_i ~ Normal(1000/r_i+z, sigma_i^2)`；所有源经同一 `z` 耦合，逐源输出是全目录联合后验的边缘，绝不是固定零点均值后的独立拟合。`sigma_mas=0` 按点质量极限处理。

数值在 `u=1000/r`、`s=ln u` 换元后用确定性自适应 Simpson（对数域相对峰值取指数、全局误差预算、尾部相对对数密度截断 -46）与零点嵌套复合 Simpson 完成；分位数由混合 CDF 二分反演，pair 概率按 `r_n<r_f ⇔ u_n>u_f` 用生存函数积分。距离分位数收敛阈值 8e-5 pc（交付容差 1e-4 pc），pair 概率误差 ≤1e-5；非有限值、非正归一化常数或触顶不收敛抛 `ValueError`，不返回部分结果。输出版本化 `star-crossmatch-distance-catalog/1`，含零点后验 `mean_mas/sigma_mas`、按 source id 排序的 `median_pc/ci_lower_pc/ci_upper_pc`（等尾 68.27%，分位点 0.15865/0.5/0.84135，有限、为正、严格有序）、按 pair id 排序的 `probability`。完整字段、数值策略与独立 oracle 交叉验证见 [整批视差贝叶斯距离契约](distance_catalog.md)。
