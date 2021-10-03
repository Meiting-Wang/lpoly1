clear all
cls

* 生成测试数据
set seed 123456
set obs 10000
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
egen x = fill(-1 `=-1+2/(50-1)')
replace x = . if _n >50 //只在50个点上进行估计

* 使用命令估计
lpoly1 Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian)

* 使用官方命令进行估计
lpoly Y X, at(x) bwidth(0.1) degree(3) kernel(gaussian) generate(beta0_cheak) nograph

* 上面两个命令的结果的对照
gen diff = abs(beta0_cheak - beta0)
list x beta0* diff in 1/20, sep(10)
sum x beta0* diff
