log using assignment1.log, replace

clear all
set seed 12345

set obs 3200
gen x = runiform()
gen z = runiform()

foreach i of numlist 40 120 400 3200 {

	tempfile mc`i'
	local j = 1
	while `j' <= 2000{
		quietly{
			preserve
			sample `i', count
			noisily di "." _continue
			
			gen e = (invnorm(uniform()))^2
			sum e
			replace e = e - r(mean)
			
			gen y = 1.2 + 1.2*x + 0.8*z + e
			
			reg y x z
			
			matrix b = e(b)
			svmat b, name(ols`i')
			collapse (mean) ols*
			
			if `j' > 1{
				append using `mc`i''
			}
			save `mc`i'', replace
			restore
			local j = `j' + 1
		}
	}
}

use `mc40', clear
append using `mc120'
append using `mc400'
append using `mc3200'
set more on

summarize

*drop x* 
*drop den*

twoway (kdensity ols401, xline(1.2)) (kdensity ols32001)
twoway (kdensity ols402, xline(0.8)) (kdensity ols32002)
twoway (kdensity ols403, xline(1.2)) (kdensity ols32003)


/* As the sample size increases, the distribution of the OLS estimates
becomes more concentrated around the true parameter values.
The empirical distributions become approximately normal despite
the skewed χ² error term.
This confirms the asymptotic normality of OLS.*/

*drop mse*
gen mse401 = (ols401 - 1.2)^2
gen mse32001 = (ols32001 - 1.2)^2
sum mse401
sum mse32001
di .6651588/.0076473

/*The MSE decreases as N increases, showing that the estimator
converges to the true parameter value.*/

*prove √N–consistent
sum ols401
scalar sd401 = r(sd)

sum ols1201 
scalar sd1201 = r(sd)

sum ols4001
scalar sd4001 = r(sd)

sum ols32001
scalar sd32001 = r(sd)

clear
set obs 4

gen N = 40 in 1
replace N = 120  in 2
replace N = 400  in 3
replace N = 3200 in 4

gen sd = .
replace sd = sd401   in 1
replace sd = sd1201  in 2
replace sd = sd4001  in 3
replace sd = sd32001 in 4

gen ln_sd = ln(sd)
gen ln_N  = ln(N)

reg ln_sd ln_N

/*The estimated slope coefficient in the regression of ln(sd) on ln(N)
is approximately -0.5.
This confirms that the standard deviation of the OLS estimator
decreases at rate 1/sqrt(N), consistent with √N-consistency.*/

*RESET
clear all
set seed 12345
set obs 500

gen x1 = runiform()
gen x2 = runiform()
gen x3 = runiform()

foreach i of numlist 30 200{
tempfile mc`i'

local j = 1
while `j' <= 4000 {
    quietly {
        preserve
        sample `i', count
        noisily di in yellow "." _continue
        
        gen e = invnorm(uniform())

        gen y = 1 + x1 + x2 + x3 + e
        
        regress y x1 x2 x3
        
        * RESET test (manual version)
        predict y_hat, xb
        gen y_hat2 = y_hat^2
        gen y_hat3 = y_hat^3
        
        regress y x1 x2 x3 y_hat2 y_hat3
        
        * F-test for joint significance
        test y_hat2 y_hat3
        
        gen p_reset = r(p)
        
        gen r_10 = p_reset < 0.10
        gen r_5  = p_reset < 0.05
        gen r_1  = p_reset < 0.01
        
        collapse (mean) r_*
        
        if `j' > 1 {
            append using `mc`i''
        }
        
        save `mc`i'', replace
        restore
        local j = `j' + 1
    }
}
}
use `mc30', clear
summarize r_10 r_5 r_1

use `mc200', clear
summarize r_10 r_5 r_1

/*For N = 30, the empirical rejection frequencies are close to
the nominal significance levels but may show small deviations
due to finite sample effects.
For N = 200, the empirical rejection rates are very close to
the theoretical levels (10%, 5%, 1%).
This indicates that the RESET test has correct asymptotic size.*/

log close