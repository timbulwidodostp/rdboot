program rdboot
		version 14.2
		set more off
		syntax varlist(min=2 numeric) [if] [in] [, c(real 0) fuzzy(varname numeric)]
		marksample touse
		tokenize `varlist'		
		tempvar x yv
		gen `x' =`2'-`c'
		gen `yv' = `1'
		
		** common parts for both sharp and fuzzy ** 
		quietly summarize `x' if `touse'
		local obs = r(N)
		set seed 98034
		mata:st_matrix("epsilon",invnormal(uniform(2500,`obs')))
		mat ulist = J(201,1,0)
		forvalues i=1/201 {
		mat ulist[`i',1]=-1+0.01*(`i'-1)
		} 
		scalar ulist_interval = 0.01

		mat Gamma_right = J(3,3,0)
		mat Gamma_left = J(3,3,0)
		mat Psi_right = J(3,3,0)
		mat Psi_left = J(3,3,0)
		mat theta12_right = J(2,1,0)
		mat theta12_left = J(2,1,0)	
				
		forvalues i=1/201 {
			scalar right = (ulist[`i',1] > 0)
			scalar left = (ulist[`i',1] < 0)
			mat r2 = (1,ulist[`i',1],ulist[`i',1]^2)
			mat ulist_i = ulist[`i',1]
			mata:K1("ulist_i")
			mat Gamma_right = Gamma_right + r2' * r2 * kernel * ulist_interval * right
			mat Gamma_left = Gamma_left + r2' * r2 * kernel * ulist_interval * left
			mat Psi_right = Psi_right + r2' * r2 * (kernel[1,1])^2 * ulist_interval * right
			mat Psi_left = Psi_left + r2' * r2 * (kernel[1,1])^2 * ulist_interval * left
			mat theta12_right = theta12_right + r2[1,1..2]' * ulist[`i',1]^2 * kernel * ulist_interval * right
			mat theta12_left = theta12_left + r2[1,1..2]' * ulist[`i',1]^2 * kernel * ulist_interval * left
		}
		
		
		
		quietly summarize `x' if `touse'
		scalar silver_h = 1.06 * r(sd) / `obs'^(1/5)
				
		
		mata:K("`x'","`touse'","silver_h")
		mata:st_numscalar("kernelsum", colsum(rowsum(st_matrix("kernel"))))
		scalar fx0 = kernelsum/`obs'/silver_h
		
		** Sharp RDD case ** 
		if "`fuzzy'" == "" {
				scalar h1 = (Gamma_right[2,2]^2 / (4*Psi_right[2,2]))^(1/5) / fx0 * `obs'^(-1/5)
				mata:estimate_right("`x'","`yv'","`touse'","h1")
				mat alphaU1right = alpharight
				
				mata:estimate_left("`x'","`yv'","`touse'","h1")
				mat alphaU1left = alphaleft
				
				mata:estimate_variance_right("`x'","`yv'","`touse'","h1","alphaU1right")
				scalar sigmaU1right = sigmaright
				
				mata:estimate_variance_left("`x'","`yv'","`touse'","h1","alphaU1left")
				scalar sigmaU1left = sigmaleft
				
				** bias and variance ** 
				mat biasU1 = inv(Gamma_right[1..2,1..2]) * theta12_right * alphaU1right[3,1]/(h1^2/2)/2
				mat biasU1 = biasU1 - inv(Gamma_left[1..2,1..2]) * theta12_left * alphaU1left[3,1]/(h1^2/2)/2
				scalar biasU1 = biasU1[1,1]
				
				mat varU1 = sigmaU1right*inv(Gamma_right[1..2,1..2])*Psi_right[1..2,1..2]*inv(Gamma_right[1..2,1..2])
				mat varU1 = varU1 + sigmaU1left*inv(Gamma_left[1..2,1..2])*Psi_left[1..2,1..2]*inv(Gamma_left[1..2,1..2])
				scalar varU1 = varU1[1,1]/fx0
				
				scalar h1 = (varU1/(biasU1^2))^0.2 * `obs'^(-0.2)
				
				scalar h1 = h1 * `obs'^(-1/20)
				
				mata:estimate_right("`x'","`yv'","`touse'","h1")
				mat alphaU1right = alpharight
				
				mata:estimate_left("`x'","`yv'","`touse'","h1")
				mat alphaU1left = alphaleft
				
				local tau1 = (alphaU1right[1,1] - alphaU1left[1,1])
				
				mata:r3(`obs',"`x'","h1","`yv'","alphaU1right","alphaU1left","`touse'","Gamma_right","fx0","Gamma_left","epsilon")
				
				tempname bootstr
				quietly svmat boot_tau1,name(`bootstr')
				quietly summarize `bootstr'
				local std = r(sd)
				local lower = `tau1'-1.96*`std'
				local upper = `tau1'+1.96*`std'
				drop `bootstr'
				local z = `tau1'/`std'
				local p = 1-normal(abs(`z'))
				
				** display results in a table ** 
				display ""
				display "Sharp RD estimates using local polynomial regression and multiplier bootstrap method"
				display ""
				di as text _col(12) " Method {c |}   Coef.    Std. Err.    z     P>|z|    [95% Conf. Interval]"
				di as text "{hline 19}{c +}{hline 60}"
				di as text _col(13) "Rdboot {c |}" as result " " %7.0g `tau1' "   "  %7.0g `std' "   "  %7.0g `z' "  " %5.0g `p' "      " %7.0g `lower' "   " %7.0g `upper'
				
				
		}
		
		
		** fuzzy RDD case ** 
		if "`fuzzy'" != "" {
				tempvar D
				gen `D' = `fuzzy'
		
					
				scalar h1 = (Gamma_right[2,2]^2 / (4*Psi_right[2,2]))^(1/5) / fx0 * `obs'^(-1/5)
				scalar b1 = (Gamma_right[2,2]^2 / (4*Psi_right[2,2]))^(1/5) / fx0 * `obs'^(-1/5)
				
				mata:estimate_right("`x'","`yv'","`touse'","h1")
				mat alphaU1right = alpharight
				
				mata:estimate_left("`x'","`yv'","`touse'","h1")
				mat alphaU1left = alphaleft
								
				mata:estimate_right("`x'","`D'","`touse'","b1")
				mat alphaT1right = alpharight
				
				mata:estimate_left("`x'","`D'","`touse'","b1")
				mat alphaT1left = alphaleft

				mata:estimate_variance_right("`x'","`yv'","`touse'","h1","alphaU1right")
				scalar sigmaU1right = sigmaright
				
				mata:estimate_variance_left("`x'","`yv'","`touse'","h1","alphaU1left")
				scalar sigmaU1left = sigmaleft
				
				mata:estimate_variance_right("`x'","`D'","`touse'","b1","alphaT1right")
				scalar sigmaT1right = sigmaright
				
				mata:estimate_variance_left("`x'","`D'","`touse'","b1","alphaT1left")
				scalar sigmaT1left = sigmaleft
		
				** bias and variance ** 
				mat biasU1 = inv(Gamma_right[1..2,1..2]) * theta12_right * alphaU1right[3,1]/(h1^2/2)/2
				mat biasU1 = biasU1 - inv(Gamma_left[1..2,1..2]) * theta12_left * alphaU1left[3,1]/(h1^2/2)/2
				scalar biasU1 = biasU1[1,1]
				
				mat varU1 = sigmaU1right*inv(Gamma_right[1..2,1..2])*Psi_right[1..2,1..2]*inv(Gamma_right[1..2,1..2])
				mat varU1 = varU1 + sigmaU1left*inv(Gamma_left[1..2,1..2])*Psi_left[1..2,1..2]*inv(Gamma_left[1..2,1..2])
				scalar varU1 = varU1[1,1]/fx0
				
				scalar h1 = (varU1/(biasU1^2))^0.2 * `obs'^(-0.2)
				
				mat biasT1 = inv(Gamma_right[1..2,1..2]) * theta12_right * alphaT1right[3,1]/(b1^2/2)/2
				mat biasT1 = biasT1 - inv(Gamma_left[1..2,1..2]) * theta12_left * alphaT1left[3,1]/(b1^2/2)/2
				scalar biasT1 = biasT1[1,1]
				mat varT1 = sigmaT1right*inv(Gamma_right[1..2,1..2])*Psi_right[1..2,1..2]*inv(Gamma_right[1..2,1..2])
				mat varT1 = varT1 + sigmaT1left*inv(Gamma_left[1..2,1..2])*Psi_left[1..2,1..2]*inv(Gamma_left[1..2,1..2])
				scalar varT1 = varT1[1,1]/fx0
				
				scalar b1 = ((1/4) * varT1/(biasT1^2))^(1/5) * `obs'^(-1/5)
				
				** generate h1 and b1 **
				
				scalar h1 = h1 * `obs'^(-1/20)
				scalar b1 = b1 * `obs'^(-1/20)
				
				mata:estimate_right("`x'","`yv'","`touse'","h1")
				mat alphaU1right = alpharight
				
				mata:estimate_left("`x'","`yv'","`touse'","h1")
				mat alphaU1left = alphaleft
									
				mata:estimate_right("`x'","`D'","`touse'","b1")
				mat alphaT1right = alpharight
					
				mata:estimate_left("`x'","`D'","`touse'","b1")
				mat alphaT1left = alphaleft

				local tau1 = (alphaU1right[1,1] - alphaU1left[1,1]) / (alphaT1right[1,1] - alphaT1left[1,1])
		
				** generate estimated process **
				mata:r2(`obs',"`x'","h1","b1","`yv'","alphaU1right","alphaU1left","`D'","`touse'","Gamma_right","fx0","Gamma_left","alphaT1right","alphaT1left","epsilon")

				tempname bootstr
				quietly svmat boot_tau2,name(`bootstr')
				quietly summarize `bootstr'
				local std = r(sd)
				local lower = `tau1'-1.96*`std'
				local upper = `tau1'+1.96*`std'
				drop `bootstr'
				local z = `tau1'/`std'
				local p = 1-normal(abs(`z'))
				
				** display results in a table ** 
				display ""
				display "Fuzzy RD estimates using local polynomial regression and multiplier bootstrap method"
				display ""
				di as text _col(12) " Method {c |}   Coef.    Std. Err.    z     P>|z|    [95% Conf. Interval]"
				di as text "{hline 19}{c +}{hline 60}"
				di as text _col(13) "Rdboot {c |}" as result " " %7.0g `tau1' "   "  %7.0g `std' "   "  %7.0g `z'  "  " %5.0g `p' "      " %7.0g `lower' "     " %7.0g `upper'
				
		
		
		}
		
		
end


version 14.2
mata:
void K (string scalar a,string scalar touse,string scalar b)
{
	real colvector X
	st_view(X=., ., a,touse)
	h = st_numscalar(b)
	k = 0.75*(J(length(X/h),1,1)-(X/h):^2):*((X/h):<1):*((X/h):>-1)
	st_matrix("kernel",k)
}
				

void K1 (string matrix u)
{	
	real colvector U
	U = st_matrix(u)
	a = 0.75*(J(length(U),1,1)-U:^2):*(U:<1):*(U:>-1)
	st_matrix("kernel",a)
}
				
function K2 (real vector u)
{
	a = 0.75*(J(length(u),1,1)-u:^2):*(u:<1):*(u:>-1)
	return(a)
}
				
				
void estimate_right(string scalar a,string scalar b,string scalar touse,string scalar c)
{
	real colvector X,Y
	st_view(X=., ., a,touse) 
	st_view(Y=., ., b,touse) 
	h = st_numscalar(c)
	delta=(X:>0)
	r2 = (J(length(X),1,1),X/h,(X/h):^2)'
	kdelta = diag((K2(X/h)):*delta)
	alphaU1right=cholsolve(r2*kdelta*r2',r2*kdelta*Y)
	st_matrix("alpharight",alphaU1right)
}
				
				
void estimate_left(string scalar a,string scalar b,string scalar touse,string scalar c)
{
	real colvector X,Y
	st_view(X=., ., a,touse)
	st_view(Y=., ., b,touse)
	h = st_numscalar(c)
	delta=(X:<0)
	r2 = (J(length(X),1,1),X/h,(X/h):^2)'
	kdelta = diag((K2(X/h)):*delta)
	alphaU1left=cholsolve(r2*kdelta*r2',r2*kdelta*Y)
	st_matrix("alphaleft",alphaU1left)
}
				
				
				
void estimate_variance_right(string scalar a,string scalar b,string scalar touse,string scalar c, string matrix alpha1)
{
	real colvector X,Y,Z
	st_view(X=., ., a,touse)
	st_view(Y=., ., b,touse)
	h = st_numscalar(c)
	Z = st_matrix(alpha1)
	delta=(X:>0)
	r2 = (J(length(X),1,1),X/h,(X/h):^2)'
	kdelta = diag(K2(X/h):*delta)
	mu = r2'*Z
	epsilon = Y - mu
	sigmaU1right = (epsilon'*kdelta*epsilon/sum(kdelta))^0.5
	st_numscalar("sigmaright",sigmaU1right)
}
				
void estimate_variance_left(string scalar a,string scalar b,string scalar touse,string scalar c, string matrix alpha1)
{
	real colvector X,Y,Z
	st_view(X=., ., a,touse)
	st_view(Y=., ., b,touse)
	h = st_numscalar(c)
	Z = st_matrix(alpha1)
	delta=(X:<0)
	r2 = (J(length(X),1,1),X/h,(X/h):^2)'
	kdelta = diag(K2(X/h):*delta)
	mu = r2'*Z
	epsilon = Y - mu
	sigmaU1left = (epsilon'*kdelta*epsilon/sum(kdelta))^0.5
	st_numscalar("sigmaleft",sigmaU1left)
}
				

void r2 (scalar N,string scalar a, string scalar b, string scalar c, string scalar d, string matrix alpha1, string matrix alpha2, string scalar e,string scalar touse,string matrix f, string scalar g,string matrix h,string matrix alpha3,string matrix alpha4,string matrix alpha5)
{
	real colvector x,y,D1,alphaU1right,alphaU1left,alphaT1right,alphaT1left,epsilon,Gamma_right,Gamma_left
	st_view(x=., ., a,touse)
	st_view(y=., ., d,touse)
	st_view(D1=., ., e,touse)
	alphaU1right = st_matrix(alpha1)
	alphaU1left = st_matrix(alpha2)
	alphaT1right = st_matrix(alpha3)
	alphaT1left = st_matrix(alpha4)
	epsilon = st_matrix(alpha5)
				
	Gamma_right = st_matrix(f)
	Gamma_left = st_matrix(h)
	h1 = st_numscalar(b)
	b1 = st_numscalar(c)
	fx0 = st_numscalar(g)
	right = (x:>0)
	left = (x:<0)
	support1 = ((x/h1):<1):*((x/h1):>-1)
	support2 = ((x/b1):<1):*((x/b1):>-1)
	r2h1 = (J(length(x),1,1),x/h1,(x/h1):^2)'
	r2b1 = (J(length(x),1,1),x/b1,(x/b1):^2)'
	E1Numerator_right = (y' - (alphaU1right' * r2h1)):*support1'
	E1Numerator_left = (y'- (alphaU1left' * r2h1)):*support1'
	E1Denominator_right = (D1' - (alphaT1right' * r2b1)):*support2'
	E1Denominator_left = (D1' - (alphaT1left' * r2b1)):*support2'
	E1Numerator_right_IFR = (cholsolve(Gamma_right,I(3))*r2h1)[1,1::N]:*E1Numerator_right:*(K2(x/h1))':*right'/(fx0*(N*h1)^0.5)
	E1Numerator_left_IFR = (cholsolve(Gamma_left,I(3))*r2h1)[1,1::N]:*E1Numerator_left:*(K2(x/h1))':*left'/(fx0*(N*h1)^0.5)
	E1Denominator_right_IFR = (cholsolve(Gamma_right,I(3))*r2b1)[1,1::N]:*E1Denominator_right:*(K2(x/b1))':*right'/(fx0*(N*b1)^0.5)
	E1Denominator_left_IFR = (cholsolve(Gamma_left,I(3))*r2b1)[1,1::N]:*E1Denominator_left:*(K2(x/b1))':*left'/(fx0*(N*b1)^0.5)

	
	nu1Numerator_right = epsilon * E1Numerator_right_IFR'
	nu1Numerator_left = epsilon * E1Numerator_left_IFR'
	nu1Denominator_right = epsilon * E1Denominator_right_IFR'
	nu1Denominator_left = epsilon * E1Denominator_left_IFR'
	B1Numerator_rootNh = (nu1Numerator_right - nu1Numerator_left)/(N*h1)^0.5
	B1Denominator_rootNh = (nu1Denominator_right - nu1Denominator_left)/(N*b1)^0.5
	boot_tau_fuzzy = ((alphaU1right[1]-alphaU1left[1])*B1Denominator_rootNh-(alphaT1right[1]-alphaT1left[1])*B1Numerator_rootNh )/(alphaT1right[1]-alphaT1left[1])^2
	st_matrix("boot_tau2",boot_tau_fuzzy)
	

}

void r3 (scalar N,string scalar a, string scalar b, string scalar d, string matrix alpha1, string matrix alpha2,string scalar touse,string matrix f, string scalar g,string matrix h,string matrix alpha5)
{
	real colvector x,y,alphaU1right,alphaU1left,Gamma_right,Gamma_left
	st_view(x=., ., a,touse)
	st_view(y=., ., d,touse)
	alphaU1right = st_matrix(alpha1)
	alphaU1left = st_matrix(alpha2) 
	epsilon = st_matrix(alpha5)
					
	Gamma_right = st_matrix(f)
	Gamma_left = st_matrix(h)
	h1 = st_numscalar(b)
	fx0 = st_numscalar(g)
	right = (x:>0)
	left = (x:<0)
	support = ((x/h1):<1):*((x/h1):>-1)
	r2h1 = (J(length(x),1,1),x/h1,(x/h1):^2)'
	E1Numerator_right = (y' - (alphaU1right' * r2h1)):*support'
	E1Numerator_left = (y'- (alphaU1left' * r2h1)):*support'
	E1Numerator_right_IFR = (cholsolve(Gamma_right,I(3))*r2h1)[1,1::N]:*E1Numerator_right:*(K2(x/h1))':*right'/(fx0*(N*h1)^0.5)
	E1Numerator_left_IFR = (cholsolve(Gamma_left,I(3))*r2h1)[1,1::N]:*E1Numerator_left:*(K2(x/h1))':*left'/(fx0*(N*h1)^0.5)
	
	nu1Numerator_right = epsilon * E1Numerator_right_IFR'
	nu1Numerator_left = epsilon * E1Numerator_left_IFR'
	boot_tau_sharp = (nu1Numerator_right - nu1Numerator_left)/(N*h1)^0.5
	st_matrix("boot_tau1",boot_tau_sharp)
	

}
	
				
end

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
				
