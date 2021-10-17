* Description: Do one-dimensional local polynomial regression.
* Author: Meiting Wang, Ph.D. Candidate, Institute for Economic and Social Research, Jinan University
* Email: wangmeiting92@gmail.com
* Created on Oct 1, 2021


*-------------------------主程序-----------------------------
program define lpoly1
version 16.0

syntax varlist(min=2 max=2 numeric) , at(varlist numeric max=1) [KERnel(string) BWidth(numlist max=1 >0) Degree(numlist max=1 >=0 integer) keep(namelist)]

* 默认值(设定 degree、keep、kernel、bwidth 的默认值)
if "`degree'" == "" {
	local degree = 3
}
if "`keep'" == "" {
	local keep "beta0"
	if `degree' >= 1 {
		forvalues i = 1/`degree' {
			local keep "`keep' beta`i'"
		}
	}
}
if "`kernel'" == "" {
	local kernel "gaussian"
}
if "`bwidth'" == "" {
	local bwidth = 0.1
}


* 错误输入识别
if ~ustrregexm("`keep'","^beta\d+(\s+beta\d+)*$") {
	dis as error "Option {opt keep} synatx error"
	error 9999
} //保证 keep 输入的是 "beta0 beta1 ..." 类似语句


* 定义临时变量
tempvar one X_std k index
if `degree' >= 1 {
	forvalues j = 1/`degree' {
		tempvar zc`j'
	}
}

* 获得因变量和自变量
gettoken Y X: varlist
if ustrtrim("`X'") == "`at'" {
	dis as error "Variables in {opt at(varname)} and {it:xvar} should not be the same."
	error 9999
} //at(varname) 和 xvar 中的变量不应该一样


* 设置结果储存矩阵
qui count if `at' != .
local at_num = r(N)
mata: BETA = J(`at_num',`degree'+1,.)


* 核心程序
qui gen `one' = 1
forvalues i = 1/`=_N' {
	if `at'[`i'] == . { //这里要求 `at' 的缺漏值都在尾端
		continue, break
	}

	preserve
	qui replace `X' = . if `Y' == . //解决当Y为缺漏值而X不为缺漏值产生的bug(当X为缺漏值为Y不为缺漏值时没有bug)
	qui gen `X_std' = (`X'-`at'[`i']) / `bwidth' //标准化 `X'
	local zc_st "`one'" //z column sentence
	if `degree' >= 1 {
		forvalues j = 1/`degree' {
			qui gen `zc`j'' = `X_std'^`j'
			local zc_st "`zc_st' `zc`j''"
		}
	}

	*- 针对不同的 kernel function 生成不同的 k 和 index(kernel function 的类型有：gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m)
	if "`kernel'" == "gaussian" {
		qui gen `k' = normalden(`X_std')
		qui gen `index' = 1
	}
	else if "`kernel'" == "epanechnikov" {
		qui gen `k' = 3/4 * (1 - `X_std'^2/5) / sqrt(5)
		qui gen `index' = (abs(`X_std') < sqrt(5))
	}
	else if "`kernel'" == "epan2" {
		qui gen `k' = 3/4 * (1 - `X_std'^2)
		qui gen `index' = (abs(`X_std') < 1)
	}
	else if "`kernel'" == "biweight" {
		qui gen `k' = 15/16 * (1 - `X_std'^2)^2
		qui gen `index' = (abs(`X_std') < 1)
	}
	else if "`kernel'" == "cosine" {
		qui gen `k' = 1 + cos(2*_pi*`X_std')
		qui gen `index' = (abs(`X_std') < 0.5)
	}
	else if "`kernel'" == "rectangle" {
		qui gen `k' = 0.5
		qui gen `index' = (abs(`X_std') < 1)
	}
	else if "`kernel'" == "triangle" {
		qui gen `k' = 1 - abs(`X_std')
		qui gen `index' = (abs(`X_std') < 1)
	}
	else if "`kernel'" == "parzen" {
		qui gen `k' = 4/3 - 8*`X_std'^2 + 8*(abs(`X_std'))^3
		qui replace `k' = 8 * (1-abs(`X_std'))^3 / 3 if (abs(`X_std')>0.5) & (abs(`X_std')<=1)
		qui gen `index' = (abs(`X_std') <= 1)
	}
	else if "`kernel'" == "gaussian_m" {
		qui gen `k' = normalden(`X_std')
		qui replace `k' = normalden(5) * (4*(6-abs(`X_std'))^5-6*(6-abs(`X_std'))^4+3*(6-abs(`X_std'))^3) if (abs(`X_std')>5) & (abs(`X_std')<=6)
		qui gen `index' = (abs(`X_std') <= 6)
	}
	else {
		dis as error "Invalid kernel function"
		error 9999
	}

	*- 得到 BETA 矩阵第 i 行的元素--针对 at[i]
	mata: BETA[`i',.] = one_dimen_lpoly("`zc_st'","`k'","`Y'","`index'")
	restore
}

*- 得到符合 derivative 的 BETA
mata: BETA = BETA * diag((factorial(0..`degree') :/ (J(1,`degree'+1,`bwidth'):^(0..`degree'))))

*- 生成 keep 中的变量
local beta_select_num = ustrwordcount("`keep'")
local beta_select_index = ustrregexra("`=usubinstr("`keep'","beta","",.)'","\s+",",")
mata: st_store((1,rows(BETA)),st_addvar("double",tokens("`keep'")),BETA[.,(`beta_select_index')+J(1,`beta_select_num',1)])

*- 输入参数展示
dis as text "Parameters:"
dis _col(6) as text "bwidth = {result:`bwidth'}"
dis _col(6) as text "degree = {result:`degree'}"
dis _col(6) as text "kernel = {result:`kernel'}"
dis _col(6) _s(2) as text "keep = {result:`keep'}"

*- 清除 Mata 中的内存
mata: mata clear
end



*-------------------子程序---------------------
*- 针对单个 `at' 的 local polynomial regression 估计程序 
version 16.0
mata:
real rowvector function one_dimen_lpoly(string scalar Z_temp, string scalar k_temp, string scalar Y_temp, string scalar index_temp) {
	real colvector beta

	st_view(Z=., ., tokens(Z_temp), index_temp)
	st_view(k=., ., k_temp, index_temp)
	st_view(Y=., ., Y_temp, index_temp)

	beta = (cholinv(cross(Z,k,Z)))*(cross(Z,k,Y)) //cross(Z,k,Y) = Z'diag(w)Y
	return(beta')
}
end
