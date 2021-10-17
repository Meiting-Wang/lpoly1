*-----------lpoly1 与命令 lpoly1_1 的比较-----------------

clear all
cls

* 测试参数
global h = 0.1
global p = 3
global kernel "gaussian" //one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m

* 生成测试数据
set seed 123456
set obs 10000
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
range x -1 1 50

* 使用命令估计
timer clear 1
timer on 1
lpoly1 Y X, at(x) bwidth($h) degree($p) kernel($kernel)
timer off 1
renames beta0-beta3, suffix(_tar)

* 使用 lpoly1_1 命令进行估计
timer clear 2
timer on 2
lpoly1_1 Y X, at(x) bwidth($h) degree($p) kernel($kernel)
timer off 2

* 上面两个命令的结果的对照
sum beta0*
sum beta1*
sum beta2*
sum beta3*

* 用时展示
timer list
