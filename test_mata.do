*！ description: do the local polynomial regression from Y on X at x
*! Meiting Wang
*! wangmeiting92@gmail.com
*! Oct 25, 2021

clear all
macro drop _all
cls


****** 参数设置
global n = 10000
global h = 0.1
global p = 2
global grid_num = 50
global kernel "gaussian" //one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m
global seed = 123456


******* 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
range x -1 1 $grid_num
order Y X x


****** 调入 Mata 函数
do one_dimen_lpoly.mata
mata: mata mosave one_dimen_lpoly(), replace


****** 使用 Mata 函数
timer clear 1
timer on 1
mata:
Y_var = st_data(.,"Y")
X_var = st_data(.,"X")
x_var = st_data(.,"x")
one_dimen_lpoly(Y_var,X_var,x_var,"$kernel",$h,$p)
end
timer off 1


****** 使用 lpoly1 函数
timer clear 2
timer on 2
lpoly1 Y X, at(x) kernel($kernel) bwidth($h) degree($p)
timer off 2


****** 用时展示
timer list
