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
