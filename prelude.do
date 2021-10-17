/************
模拟模型：
do the local polynomial regression of Y on X
*************/


*--------------------------准备工作----------------------
****** 参数设置
clear all
macro drop _all
cls
global n = 10000
global h = 0.1 //bandwidth
global degree = 4
global kernel "gaussian"
global seed = 123456

/*
可供选择的核函数有: gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen, gaussian_m
*/



****** 生成测试数据
if "$seed" != "" {
	set seed $seed
}
set obs $n
gen X = rnormal()
gen Y = 1 + 5*X + rnormal()
range x -1 1 50
order Y X x
save testdata.dta, replace


****** 调入 Mata 函数(这个函数能计算 one dimensional 的 lpoly，并把结果储存在 Stata 数据集的对应变量中)
mata:
void function one_dimen_lpoly(
	string scalar Y_var, //dependent variable
	string scalar X_var, //independent variable
	string scalar x_var, //grid point(the missing value must be at the end, if any.)
	string scalar kernel, //kernel function(the function written must be one of the series of functions specified in the program)
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


****** 使用函数
use testdata,clear
timer clear 1
timer on 1
mata: one_dimen_lpoly("Y","X","x","$kernel",0.1,$degree)
timer off 1

timer clear 2
timer on 2
lpoly Y X, kernel($kernel) bwidth(0.1) degree($degree) at(x) generate(beta0_cheak) nograph
// lpoly1_1 Y X, at(x) bwidth(0.1) degree($degree) kernel($kernel)
timer off 2


****** 用时展示
timer list


****** 结果对比
keep x beta0*
sum
