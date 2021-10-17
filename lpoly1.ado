* Description: Do one-dimensional local polynomial regression.
* Author: Meiting Wang, Ph.D. Candidate, Institute for Economic and Social Research, Jinan University
* Email: wangmeiting92@gmail.com
* Created on Oct 1, 2021
* Updated on Oct 17, 2021


*-------------------------主程序-----------------------------
program define lpoly1
version 16.0

syntax varlist(min=2 max=2 numeric) , at(varlist numeric max=1) [KERnel(string) BWidth(numlist max=1 >0) Degree(numlist max=1 >=0 integer) keep(namelist)]

* 默认值(设定 degree keep kernel bwidth 的默认值)
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


* 获得因变量, 自变量和格点变量
tokenize `varlist'
local Y "`1'"
local X "`2'"
local x "`at'"
if "`x'" == "`X'" {
	dis as error "Variables in {opt at(varname)} and {it:xvar} should not be the same."
	error 9999
} //以保证 at(varname) 和 xvar 中的变量不一样


* 使用子程序中的 one_dimen_lpoly() 函数，计算 one dimensional 的 lpoly, 并保存所设定的变量
qui ds
local varlist_old "`r(varlist)'"
mata: one_dimen_lpoly("Y","X","x","`kernel'",`bwidth',`degree')
keep `varlist_old' `keep'


* 输入参数展示
dis as text "Parameters:"
dis _col(6) as text "bwidth = {result:`bwidth'}"
dis _col(6) as text "degree = {result:`degree'}"
dis _col(6) as text "kernel = {result:`kernel'}"
dis _col(6) _s(2) as text "keep = {result:`keep'}"

end



*--------------------------子程序-----------------------------
/* 该 Mata 函数能计算 one dimensional 的 lpoly，并把结果储存在 Stata 数据集的对应变量中(beta0 beta1 ... betap) */
version 16.0
mata:
void function one_dimen_lpoly(
	string scalar Y_var, //dependent variable
	string scalar X_var, //independent variable
	string scalar x_var, //grid point(the missing value must be at the end, if any.)
	string scalar kernel, //kernel function(the function written must be one of gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m)
	real scalar h, //bandwidth(bandwidth needs to be greater than 0)
	real scalar p) //degree(degree needs to be a non-negative integer)
{
	/* 错误信息判断 */
	if (h<=0) {
		printf("{error:Bandwidth needs to be greater than 0}\n")
		exit(9999)
	} //保证 h>0
	if ( (p<0) | (mod(p,1)) ) {
		printf("{error:Degree should be a non-negative integer}\n")
		exit(9999)
	} //保证 p 为非负整数

	/* 得到没有缺漏值的 Y_raw X_raw x */
	Y_raw = st_data(.,Y_var)
	X_raw = st_data(.,X_var)
	x = st_data(.,x_var)
	Y_X_nomiss = (Y_raw:!=.):&(X_raw:!=.)
	Y_raw = select(Y_raw,Y_X_nomiss)
	X_raw = select(X_raw,Y_X_nomiss)
	x = select(x,x:!=.)
	
	/* 得到 q */
	q = rows(x)

	/* 计算BETA */
	BETA = J(q,p+1,.) //设置结果储存矩阵
	for (i=1; i<=q; i++) {
		X_std_raw = (X_raw:-x[i,1])/h

		/* 依据不同的核函数对 k_index X_std k 进行不同的处理(核函数的类别有: gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m) */
		if (kernel=="gaussian") {
			k_index = J(rows(X_std_raw),1,1)
			X_std = select(X_std_raw,k_index) //已经去除了缺漏值和考虑了 k_index 的 X_std(所有后面的操作都是以此为基础)
			k = normalden(X_std)
		}
		else if (kernel=="epanechnikov") {
			k_index = abs(X_std_raw):<sqrt(5)
			X_std = select(X_std_raw,k_index)
			k = 3/4 * (1 :- X_std:^2/5) / sqrt(5)
		}
		else if (kernel=="epan2") {
			k_index = abs(X_std_raw):<1
			X_std = select(X_std_raw,k_index)
			k = 3/4 * (1 :- X_std:^2)
		}
		else if (kernel=="biweight") {
			k_index = abs(X_std_raw):<1
			X_std = select(X_std_raw,k_index)
			k = 15/16 * (1 :- X_std:^2):^2
		}
		else if (kernel=="cosine") {
			k_index = abs(X_std_raw):<0.5
			X_std = select(X_std_raw,k_index)
			k = 1 :+ cos(2*pi()*X_std)
		}
		else if (kernel=="rectangle") {
			k_index = abs(X_std_raw):<1
			X_std = select(X_std_raw,k_index)
			k = J(rows(X_std),1,0.5)
		}
		else if (kernel=="triangle") {
			k_index = abs(X_std_raw):<1
			X_std = select(X_std_raw,k_index)
			k = 1 :- abs(X_std)
		}
		else if (kernel=="parzen") {
			k_index = abs(X_std_raw):<=1
			X_std = select(X_std_raw,k_index)
			k = J(rows(X_std),1,.)
			cond1 = selectindex(abs(X_std):<=0.5)
			cond2 = selectindex(abs(X_std):>0.5)
			k[cond1] = 4/3 :- 8*X_std[cond1]:^2 + 8*abs(X_std[cond1]):^3
			k[cond2] = 8 * (1:-abs(X_std[cond2])):^3 / 3
		}
		else if (kernel=="gaussian_m") {
			k_index = abs(X_std_raw):<=6
			X_std = select(X_std_raw,k_index)
			k = J(rows(X_std),1,.)
			cond1 = selectindex(abs(X_std):<=5)
			cond2 = selectindex(abs(X_std):>5)
			k[cond1] = normalden(X_std[cond1])
			k[cond2] = normalden(5) * (4*(6:-abs(X_std[cond2])):^5:-6*(6:-abs(X_std[cond2])):^4:+3*(6:-abs(X_std[cond2])):^3)
		}
		else {
			printf("{error:Wrong kernel function}\n")
			exit(9999)
		}

		Z = mm_expand(X_std,1,p+1):^(0..p)
		Y = select(Y_raw,k_index)
		beta = (cholinv(cross(Z,k,Z)))*(cross(Z,k,Y)) //cross(Z,k,Y) = Z'diag(w)Y
		BETA[i,.] = beta'
	}
	BETA = BETA * diag((factorial(0..p) :/ (J(1,p+1,h):^(0..p)))) //得到符合 derivative 的 BETA

	/* 将 BETA 转化为 Stata 内存中的变量(beta0 beta1 ... betap) */
	st_store((1,rows(BETA)), st_addvar("double",J(1,p+1,"beta"):+strofreal(0..p)), BETA)
}
end
