# Monte Carlo Study of OLS Asymptotics and RESET Test
This project is an assignment for my Advanced Applied Econometrics course at KU Leuven, which uses Monte Carlo simulations. 

## 💻 Implementation and Purpose
The implementation was carried out in Stata, with several objectives. The purposes were:
- Show that the OLS error term does not need to be normally distributed.
- Show that OLS is √N–consistent
- Use Monte Carlo simulations to study the empirical performance of the test and verify whether the Type I error rate corresponds to the theoretical significance level.
## 📚 What I Learned
The Ordinary Least Squares estimator has the Gauss–Markov theorem property of being BLUE (Best Linear Unbiased Estimator). A Monte Carlo simulation can be used to verify the consistency property. Furthermore, the sample size directly affects the standard error, which decreases at a rate of 1/√N.

This relationship can be shown as: 
ln(σ^)=α ln(N) + c

where the theoretical prediction is:
α = -1/2

Finally, the implementation of the Ramsey RESET test is useful for detecting general functional form misspecification. The RESET test verifies the hypothesis that the true functional relationship is linear by adding polynomial terms of the OLS fitted values to the regression. If these additional terms are statistically significant, this indicates model misspecification.
## 🌟 Conclusion
The task has been very useful for understanding the properties of the Ordinary Least Squares estimator and how they can be demonstrated through different methods, such as Monte Carlo simulation.
