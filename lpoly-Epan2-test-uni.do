/*
- 仅计算一个点（元素）的局部多项式回归，用以与命令 lpoly 计算的结果进行对比
- 以 Epan2 核函数为例
*/

*---------------------自己编程计算结果---------------------------
clear all
macro drop _all
cls

* 开始计时
timer clear 1
timer on 1


* 参数设置
global x = 0.1
global h = 0.1
global degree = 3

* 创建数据
set obs 10000
set seed 123456
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
gen X_std = (X-$x)/$h
save temp.dta, replace


* 生成设计矩阵 Z 所需变量 one zc1 zc2 zc3 ... 和生成核函数矩阵 K 所需的变量 k (以 Epan2 核函数为例)
use temp.dta, clear
gen one = 1
local zc_st "one"
if $degree >= 1 {
	forvalues i = 1/$degree {
		gen zc`i' = X_std^`i'
		local zc_st "`zc_st' zc`i'"
	}
}

gen k = 3/4 * (1 - X_std^2)
gen d = (abs(X_std) < 1)


* 使用 mata 编程实现 local polynomial regression
mata:
st_view(Z=., ., tokens("`zc_st'"), "d")
st_view(k=., ., "k", "d")
st_view(Y=., ., "Y", "d")

beta = (cholinv(cross(Z,k,Z)))*(cross(Z,k,Y)) //cross(Z,k,Y) = Z'diag(w)Y
st_numscalar("Yhat1",beta[1,1])
end


*--------------------使用lpoly命令计算结果---------------------------
use temp.dta, clear
gen X_grid = .
replace X_grid = $x in 1
lpoly Y X, kernel(epan2) bwidth($h) degree($degree) at(X_grid) generate(beta) nograph
scalar Yhat2 = beta[1]


*---------------------两者相互比较------------------------
scalar diff = Yhat2 - Yhat1
scalar list

* 结束计时
timer off 1
timer list 1

