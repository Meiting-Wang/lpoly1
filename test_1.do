clear all
cls

* 生成测试数据
set seed 123456
set obs 100000
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
range x -1 1 50

* 使用命令估计
timer clear 1
timer on 1
lpoly1_1 Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian)
timer off 1

* 使用官方命令进行估计
timer clear 2
timer on 2
lpoly Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian) generate(beta0_cheak) nograph
timer off 2

* 上面两个命令的结果的对照
gen diff_beta0 = abs(beta0_cheak - beta0)
sum diff_beta0

* 用时展示
timer list
