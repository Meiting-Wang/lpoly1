/************
模拟模型：
do the local polynomial regression of Y on X
*************/


*--------------------------准备工作----------------------
* 参数设置
clear all
macro drop _all
cls
global n = 10000
global h = 0.1 //bandwidth
global degree = 3
global kernel "epan2"
global seed = 123456



* 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
egen x = fill(-1 `=-1+2/(50-1)')
replace x = . if _n >50 //只在50个点上进行估计
save testdata.dta, replace


* 设置结果储存矩阵，调入 Mata 函数
mata:
BETA = J(50,$degree+1,.)
real rowvector function one_lpoly(string scalar Z_temp, string scalar k_temp, string scalar Y_temp, string scalar index_temp) {
	real colvector beta

	st_view(Z=., ., tokens(Z_temp), index_temp)
	st_view(k=., ., k_temp, index_temp)
	st_view(Y=., ., Y_temp, index_temp)

	beta = (cholinv(cross(Z,k,Z)))*(cross(Z,k,Y)) //cross(Z,k,Y) = Z'diag(w)Y
	return(beta')
}
mata describe
end



* 主程序
cls
use testdata.dta, clear
gen one = 1

forvalues i = 1/`=_N' {
	if x[`i'] == . { //这里要求 x 的缺漏值都在尾端
		continue, break
	}

	preserve
	gen X_std = (X-x[`i']) / $h //标准化 X
	local zc_st "one" //z column sentence
	if $degree >= 1 {
		forvalues j = 1/$degree {
			gen zc`j' = X_std^`j'
			local zc_st "`zc_st' zc`j'"
		}
	}
	gen k = 3/4 * (1 - X_std^2)
	gen index = (abs(X_std) < 1)

	mata: BETA[`i',.] = one_lpoly("`zc_st'","k","Y","index")
	restore
	dis in y _c "`i' "
}


mata:
st_store((1,rows(BETA)),st_addvar("double",tokens("beta0 beta1 beta2 beta3")),BETA)
end


* lpoly 命令的检验
lpoly Y X, kernel(epan2) bwidth($h) degree($degree) at(x) generate(beta0_cheak) nograph

list x beta0* in 1/50
