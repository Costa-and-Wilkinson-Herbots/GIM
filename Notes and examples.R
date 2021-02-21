### R CODE TO SUPPLEMENT COSTA & WILKINSON-HERBOTS (2021):
### "Inference of gene flow in the process of speciation: 
### efficient maximum-likelihood implementation of a generalised isolation-with-migration model".

### Before you use the code in this repository, you need to:
### - install and load the packages 'stats' (contains the function 'nlminb')
###   and 'numDeriv' (contains the function 'hessian'); 
### - if you wish to fit models to data: run the R code in the file "Functions to fit models";
### - if you wish to simulate data: run the R code in the file "Functions to simulate data".

### 1) To simulate data from the GIM model or simpler models nested within it:
#
#   Example:
#   to simulate data from scenario (i) in the paper, with mutation rate heterogeneity, using the parameters in Figure 2 and eq.(12) of the paper:

data<-sim_r.GIM(theta0=3,theta1=2,theta2=4,theta1_prime=3,theta2_prime=6,T1=4,V=4,M1_star=0.2,M2_star=0.4,M1_prime_star=0.04,M2_prime_star=0.08,nums=c(10000,10000,20000),alpha=10,lambda=10)

#   Result:
#   data$x1 will contain 10000 observations of the number of nucleotide differences between two sequences from species 1;
#   data$x2 will contain 10000 observations of the number of nucleotide differences between two sequences from species 2;
#   data$x3 will contain 20000 observations of the number of nucleotide differences between a sequence from species 1 and a sequence from species 2;
#   data$r1, data$r2 and data$r3 will contain the relative mutation rates (or mutation scalars) for the loci in data$x1, data$x2 and data$x3, respectively;
#   data$t1, data$t2 and data$t3 will contain the corresponding coalescence times;
#   data$n1, data$n2 and data$n3 will contain the numbers of observations generated (10000, 10000, 20000 in this example).
# 
#   Notes:
#   To simulate without mutation rate heterogeneity, use sim.GIM instead.
#   To simulate using the original model parameters (see Figure 1d), use sim_r.GIM.original (with mutation rate heterogeneity) 
#   or sim.GIM.original (without mutation rate heterogeneity) instead.
#
#   For further information, please see the notes in the file "Functions to simulate data".

### 2) To fit a model to data:
#
#   Example:
#   to fit the GIM model to the data simulated above:

fit.GIM<-gim(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"GIM")

#   Result:
#   fit.GIM$estimates contains the vector of ML estimates (easier to use in further computations than the labelled version below);
#   fit.GIM$estimates_labelled contains a list containing the ML estimate of each parameter, labelled by its name for clarity;
#   fit.GIM$max_lnL contains the maximized loglikelihood;
#   fit.GIM$AIC contains the AIC score of the fitted model;
#   fit.GIM$fit contains the result of nlminb, including information regarding convergence.
#   The function also prints the ML estimates, maximized loglikelihood and AIC score.
#
#   Note: the parameters are as defined by the reparameterisation in eq.(12) of the paper (see also Figure 2).

#   to fit the IIM model to the data simulated above:

fit.IIM<-gim(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"IIM")

#   to fit the secondary contact model to the data simulated above:

fit.SEC<-gim(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"SEC")

#   to fit the isolation model to the data simulated above:

fit.ISO<-gim(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"ISO")

#   Note:
#   Default starting values and lower bounds for nlminb are specified in the file "Functions to fit models".
#   Depending on the data, these may need to be changed to more suitable values.
#   Example:
#   to fit the GIM model to the data simulated above, specifying a starting value (iv) and lower bound (lb):

fit.GIM<-gim(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"GIM",iv=c(rep(1,5),rep(2,2),rep(0.5,4)),lb=c(rep(0.05,5),rep(0,6)))

#   Further notes:
#    - The loglikelihood may have multiple local maxima, so ideally a number of different starting values should be tried.
#    - When fitting an 'incorrect' model (for example, fitting a model of secondary contact to data generated from a model with decreasing gene flow), 
#      the likelihood maximization function may try smaller and smaller values of the population size parameters, resulting in infinite or NaN values; 
#      if this happens, we suggest increasing the lower bounds on the population size parameters in the optional input vector lb; 
#      trying different starting values may also help.
#    - The estimates are meaningful only if 'nlminb' returns "relative convergence" or "absolute convergence". If not, other starting values should be tried.

### 3) To perform, for example, a Likelihood Ratio test of the isolation model (H0) vs. the IIM model (H1) for the data simulated above, 
#   using the mixture \frac{1}{2} \chi^2_1 + \frac{1}{2} \chi^2_2 as the null distribution, compute the LRT statistic and the p-value:

LRT <- 2*(fit.IIM$max_lnL - fit.ISO$max_lnL)
p <- (pchisq(LRT,df=1,lower.tail=FALSE)+pchisq(LRT,df=2,lower.tail=FALSE))/2

#   Alternatively, compute the appropriate quantile of the null distribution. 
#   For example, to perform the LRT test above at a significance level of 2.5%:

q <- quantile(c(rchisq(10**8,df=1),rchisq(10**8,df=2)),probs=0.975)

### 4) To compute Wald standard errors and the covariance matrix of the ML estimates:
#   
#   Example: to compute standard errors of the ML estimates of the parameters of the GIM model for the data simulated above:

fit.GIM.se<-gim_se(x1=data$x1,x2=data$x2,x3=data$x3,r1=data$r1,r2=data$r2,r3=data$r3,"GIM",fit.GIM$estimates)

#   Result:
#   fit.GIM.se$se contains the vector of standard errors of the ML estimates;
#   fit.GIM.se$Cov contains the covariance matrix of the ML estimates.
#   The function also prints the standard errors of the ML estimates.
# 
#   Note: 
#   Our code to compute the Wald standard errors and covariance matrix will not work for a given model if one or more parameters are estimated to be zero. 
#   An alternative approach in such cases is to reduce the model by fixing at zero the parameter(s) which have been estimated to be zero, and to compute 
#   the Wald standard errors for the parameters of the reduced model. As the data will have the same likelihood under both models and the reduced model has 
#   fewer parameters, the latter model is preferable.

#   For further information, please see the notes in the file "Functions to fit models".



