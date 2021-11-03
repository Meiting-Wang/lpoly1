clear all
macro drop _all
cls


* 参数设置
global n = 10000
global h = 0.1
global p = 2
global grid_num = 50
global kernel "gaussian" //one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m
global seed = 123456


* 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
range x -1 1 $grid_num
order Y X x


* 使用命令估计
timer clear 1
timer on 1
lpoly1 Y X, at(x) kernel($kernel) bwidth($h) degree($p)
timer off 1


* 使用官方命令进行估计
timer clear 2
timer on 2
lpoly Y X, at(x) kernel($kernel) bwidth($h) degree($p) generate(beta0_cheak) nograph
timer off 2


* 上面两个命令的结果的对照
sum beta0*


* 用时展示
timer list
