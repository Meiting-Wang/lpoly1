mata:
real matrix function one_dimen_lpoly(
	real colvector Y_var, //dependent variable
	real colvector X_var, //independent variable
	real colvector x_var, //grid point(the missing value must be at the end, if any.)
	string scalar kernel, //kernel function(the function written must be one of gaussian_m, gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen)
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

	/* 得到没有缺漏值的 Y_raw X_raw x, (Y_raw X_raw) 匹配 */
	Y_X_nomiss = (Y_var:!=.):&(X_var:!=.)
	Y_raw = select(Y_var,Y_X_nomiss)
	X_raw = select(X_var,Y_X_nomiss)
	x = select(x_var,x_var:!=.)
	
	/* 得到 q */
	q = rows(x)

	/* 计算BETA */
	BETA = J(q,p+1,.) //设置结果储存矩阵
	for (i=1; i<=q; i++) {
		X_std_raw = (X_raw:-x[i,1])/h

		/* 依据不同的核函数对 k_index X_std k 进行不同的处理(核函数的类别有: gaussian_m, gaussian, epanechnikov, epan2, biweight, cosine, rectangle, triangle, parzen) */
		if (kernel=="gaussian") {
			k_index = J(rows(X_std_raw),1,1)
			X_std = select(X_std_raw,k_index) //已经去除了缺漏值和考虑了 k_index 的 X_std, 所有后面的操作都是以此为基础(下同)
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

	/* 返回 BETA */
	return(BETA)
}
end
