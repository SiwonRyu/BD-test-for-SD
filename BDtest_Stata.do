
/* Stochastic Dominance Test of Barrett-Donald(2003) */

cap program drop bdtest
program bdtest ,eclass 
	syntax varlist(min=2 max=2) [,b(integer 10) grid(integer 10) s(integer 10)]
	ereturn clear
	mata : bdtest("`varlist'", B = `b', grid_num = `grid', s=`s')  
	ereturn scalar test_stat = bd
	ereturn matrix pvalue = pv
end
qui{
mata
mata clear
real matrix bdtest(varlist, B, grid_num, s){	
	st_view(data=.,.,varlist)
	x1 = data[,1]
	x2 = data[,2]
	n1 = sum(rowmissing(data[,1]):==0)
	n2 = sum(rowmissing(data[,2]):==0)
	x_grid = rangen(min(x1\x2),max(x1\x2),grid_num)

	bd = stat(x1,x2,x_grid,s)
	pval_m1 = mean(multiplier1(x1,x2,x_grid,s,B) :> J(B,1,bd))
	pval_m2 = mean(multiplier2(x1,x2,x_grid,s,B) :> J(B,1,bd))
	pval_b1 = mean(boot1(x1,x2,x_grid,s,B) :> J(B,1,bd))
	pval_b2 = mean(boot1(x1,x2,x_grid,s,B) :> J(B,1,bd))
	pval_b3 = mean(boot3(x1,x2,x_grid,s,B) :> J(B,1,bd))
	pval = (pval_m1 \ pval_m2 \ pval_b1 \ pval_b2 \ pval_b3)
	
	st_numscalar("bd", bd)
	st_matrix("pv", pval)
		
	/* Print Result */
	printf("{txt} Barrett-Donald Test for Stochastic Dominance \n")
	printf("{txt} * H0 : x1 FSD x2 \n")
	printf("{txt} * Number of grid \t= %5.0g \n",grid_num)
	printf("{txt} * Number of bootstrap \t= %5.0g\n",B)
	printf("{txt} * Dominance order \t= %5.0g\n\n",s)
	
	printf("{txt}{space 25}{c |} {space 3} BD {space 5} p-value {space 5}\t\n")
	printf("{hline 25}{c +}{hline 26}\t\n")
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Multiplier method1",bd,pval_m1)
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Multiplier method2",bd,pval_m2)
	printf("{hline 25}{c +}{hline 26}\t\n")
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Bootstrap method1",bd,pval_b1)
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Bootstrap method2",bd,pval_b2)
	printf("{txt} %20s \t {c |} {space 2} %5.4f {space 3} %5.4f \t\n","Bootstrap method3",bd,pval_b3)	
	
}
real matrix oper(X,z,s){
	Z = J(rows(X),1,z')
	op = J(rows(X),rows(z),0)
		
	for (j=1; j<=rows(z); j++){
		op[,j] = (X :<= Z[,j]):*(Z[,j]-X):^(s-1):/ factorial(s-1)
	}
	return(op)
}
real matrix ecdf(X,z,s){
	return(mean(oper(X,z,s)))
}
real matrix stat(x1,x2,z,s){
	n1 = rows(x1)
	n2 = rows(x2)
	stat = sqrt(n1*n2/(n1+n2))*max(ecdf(x1,z,s)-ecdf(x2,z,s))
	return(stat)
}
// Multiplier Methods to calculate p-value
real matrix opJ(X,z,s){
	D = oper(X,z,s)-J(rows(X),1,ecdf(X,z,s))
	U = J(1,cols(z),rnormal(rows(X),1,0,1))
	jv = sqrt(rows(X))*mean(D :* U)
	return(jv)
}
real multiplier1(x1,x2,z,s,B){
	n1 = rows(x1)
	n2 = rows(x2)
	lambda = n2/(n1+n2)
	m1 = J(B,1,0)
	for(j=1;j<=B;j++){
	m1[j,] = max(sqrt(lambda)*opJ(x1,z,s) - sqrt(1-lambda)*opJ(x2,z,s))
	}
	return(m1)
}
real multiplier2(x1,x2,z,s,B){
	n1 = rows(x1)
	n2 = rows(x2)
	m2 = J(B,1,0)
	for(j=1;j<=B;j++){
	m2[j,] = max(opJ(x2,z,s))
	}
	return(m2)
}
// Bootstrap Methods to calculate p-value
real boot1(x1,x2,z,s,B){
	n2 = rows(x2)
	b1 = J(B,1,0)
	for(j=1;j<=B;j++){
	index  = floor((n2-1):*runiform(n2,1):+1)
	b1[j,] = sqrt(n2)*max( ecdf(x2[index],z,s)-ecdf(x2,z,s) )
	}
	return(b1)
}
real boot2(x1,x2,z,s,B){
	n1 = rows(x1)
	n2 = rows(x2)
	n = n1+n2
	x = (x1\x2)
	b2 = J(B,1,0)
	for(j=1;j<=B;j++){
	index  = floor((n-1):*runiform(n,1):+1)
	xp = x[index]
	x1p = xp[1..n1,1]
	x2p = xp[n1+1..n2,1]
	b2[j,] = sqrt(n2*n1/(n1+n2))*max( ecdf(x1p,z,s)-ecdf(x2p,z,s) )
	}
	return(b2)
}
real boot3(x1,x2,z,s,B){
	n1 = rows(x1)
	n2 = rows(x2)
	b3 = J(B,1,0)
	for(j=1;j<=B;j++){
	index1  = floor((n1-1):*runiform(n1,1):+1)
	index2  = floor((n2-1):*runiform(n2,1):+1)
	b3[j,] = sqrt(n2*n1/(n1+n2))*max( ecdf(x1[index1],z,s)-ecdf(x1,z,s)-(ecdf(x2[index2],z,s)-ecdf(x2,z,s)) )
	}
	return(b3)
}

end
}

/******************************************************************************/
* Test
* H0 : x1 FSD x2 // H1 : x1 not FSD x2
* x1 ~ N(-0.5,1), x2 ~ N(0,1) : x2 FSD,SSD x1
* True : reject H0 for s = 1,2

import excel "C:\data_temp\ex.xlsx",sheet("Sheet1") first clear
// clear
// global n 100
// set obs $n
//
// gen x1 = rnormal(-0.5,1)
// gen x2 = rnormal(0,1)

bdtest x1 x4,b(500) grid(50) s(1)
*ereturn list



/* Simulation
capture program drop bdsim
program bdsim,rclass
clear
global n 100
set obs $n

gen x1 = rnormal(-0.5,1)
gen x2 = rnormal(0,1)

bdtest x1 x2,b(100) grid(50) s(1)
mat pv = e(pvalue)
return scalar pvm1 = pv[1,1] < 0.05
return scalar pvm2 = pv[2,1] < 0.05
return scalar pvb1 = pv[3,1] < 0.05
return scalar pvb2 = pv[4,1] < 0.05
return scalar pvb3 = pv[5,1] < 0.05
end

simulate pvm1 = r(pvm1) pvm2 = r(pvm2) pvb1 = r(pvb1) pvb2 = r(pvb2) pvb3 = r(pvb3), reps(100) seed(123447) : bdsim
su
