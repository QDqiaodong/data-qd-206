# 天文目录匹配契约

使用者是天文观测目录整理人员，输入reference（1..1000项id/ra/dec/magnitude）与observations（1..2000项id/ra/dec/sigma），threshold为正数。各数组id唯一ASCII；数值是最多9位整数和6位小数的十进制字符串。ra取[0,360)度，dec取[-90,90]度，sigma为正的各向同性角误差（度）。magnitude只保留为来源元数据，不参与筛选。

角距离按球面haversine计算，采用Python双精度三角函数；可跨赤经0/360与极区。距离<=threshold*sigma才是候选，边界按实际浮点结果比较，无额外epsilon。distance为浮点值字符串，不声称精确有理角距离。候选按距离、reference id排序；恰一候选才selected，多个候选ambiguous且selected=null，零候选unmatched。各观测独立，一个参考源可被多个观测命中，不进行全局一对一匹配。

这不是误差椭圆、相机重建或天文图像叠加；没有自行运动/历元换算、球面协方差、亮度联合推断。后续历元传播和源关联可接入reference/observation身份与候选证据，不用改造为通用报表。样例仅为合成角坐标。process输入错误抛ValueError；io.load严格拒绝重复JSON键，io.publish原子新建目标，已有文件或链接拒绝覆盖。
