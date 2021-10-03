cls
mata
I = I(2)
I
end


mata: I = I(3)
mata: I

sysuse auto.dta, clear
sum price

mata
stata("sum price")
end

h mata m4 matrix

mata
A = (1,2,3 \ 4,5,6)
A
e = e(2,5)
e
J = J(2,3,999)
J
end

mata: a = 8::15
mata: a'


set seed 123456
mat B = matuniform(2,3)
mat list B

mata
C = st_matrix("B")
C

mata_A = (4,5,6\7,8,9)
mata_A
st_matrix("stata_A",mata_A)
end
mat list stata_A

mata
A = (1,2,3\4,5,6)
B = (4,7,5\9,5,1)
C = (2,3\7,4\8,5)
A
B
C

sqrt(B)
end

mata
A = (1,5,4\7,8,9\1,4,2)
A
A[|2,1 \ 3,3|]
end


sysuse auto.dta, clear
sum price
dis r(sum)
mata
mata clear
void calcsum(varname, resultissum)
{
	st_view(x=., ., varname)
	resultissum = colsum(x)
}
sum = .
calcsum("price", sum)
sum
end


sysuse auto.dta, clear
sum price
dis r(sum)

mata
mata clear

void calcsum(varname)
{
	st_view(x=., ., varname)
	st_numscalar("summ",colsum(x))
}

calcsum("price")
end
ret list
scalar list





mata
A = 1.123456
st_numscalar("xx",A)

end
scalar list
ret list



clear
input x1 x2 x3 index
1 3 4 0
2 5 1 1
5 1 2 0
5 4 1 1
end


mata
st_view(Z=.,.,tokens("x1 x2 x3"))
st_view(d=.,.,tokens("index"))
Z
d

selectindex(d)
xx = select(Z,d)
xx
isview(xx)

st_select(A, Z, d)
A
isview(A)
end


*-----------------------情形一-----------------------------
clear
input x1 x2 x3 index
1 3 4 0
2 5 1 1
5 1 2 0
5 4 1 1
end


mata
st_view(d=.,.,tokens("index"))
d = selectindex(d)
st_view(Z=.,d,tokens("x1 x2 x3"))
d
Z
isview(Z)
isview(d)
end


*--------------------------情形二-----------------------------
clear
input x1 x2 x3 index
1 3 4 0
2 5 1 1
5 1 2 0
5 4 1 1
end

mata
st_view(Z=.,.,tokens("x1 x2 x3"))
st_view(d=.,.,tokens("index"))
st_select(Z_target=., Z, d)
Z
Z_target
isview(Z)
isview(d)
isview(Z_target)

st_select(XX=.,"x1",d)
end




clear
input x1 x2 x3 index
1 3 4 0
2 5 1 1
5 1 2 0
5 4 1 1
end
gen xx = x1^3


clear
input x1 x2 x3 index
1 3 4 0
2 5 1 1
5 1 2 0
5 4 1 1
end
l
levelsof x3
ret list
