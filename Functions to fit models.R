### FUNCTIONS TO FIT THE MODELS IN COSTA & WILKINSON-HERBOTS (2021): 
### "Inference of gene flow in the process of speciation: efficient maximum-likelihood implementation of a generalised isolation-with-migration model".

### This file contains the functions needed to fit the following models:
### - GIM model ("GIM");
### - IIM model ("IIM");
### - secondary contact model ("SEC");
### - isolation model ("ISO").
### Data should consist of the numbers of nucleotide differences between one pair of DNA sequences from each of a large number of independent loci.


### Function to fit a model to the data:
#
# Required input:
# x1 = vector containing the numbers of nucleotide differences between pairs of sequences from species 1
#     (this vector should contain one entry for each locus at which two sequences from species 1 are compared);
# x2 = vector containing the numbers of nucleotide differences between pairs of sequences from species 2
#     (this vector should contain one entry for each locus at which two sequences from species 2 are compared);
# x3 = vector containing the numbers of nucleotide differences between pairs of sequences from different species
#     (this vector should contain one entry for each locus at which a sequence from species 1 and a sequence from species 2 are compared);
# r1 = vector containing the relative mutation rates of the loci in x1
#     (this vector must have the same length as x1 and all its entries must be non-zero);
# r2 = vector containing the relative mutation rates of the loci in x2
#     (this vector must have the same length as x2 and all its entries must be non-zero);
# r3 = vector containing the relative mutation rates of the loci in x3
#     (this vector must have the same length as x3 and all its entries must be non-zero);
# model: "GIM", "IIM", "SEC" or "ISO";
# 
# Optional input (if not specified, default values will be used):
# iv = starting value (vector of parameters) for the likelihood maximization function;
# lb = vector of lower bounds for the parameters
#     (for the population size parameters, a lower bound > 0 must be provided;
#      for the time and migration rate parameters, the lower bound is normally set to zero.)
#
# The function returns the following:
# fit: the result of nlminb, including information regarding convergence;
# estimates = the vector of ML estimates (easier to use in further computations than the labelled version below);
# estimates_labelled = a list containing the ML estimate of each parameter, labelled by its name (for clarity);
# max_lnL = the maximized loglikelihood;
# AIC = the AIC score of the fitted model.
#
# The function also prints the ML estimates, maximized loglikelihood and AIC score.
# 
# Notes:
# (i) the parameters are as defined by the reparameterisation in eq.(12) of the paper (see also Figure 2);
# (ii) if the 'relative mutation rates' entered don't average to 1, this is not a problem; they will be automatically rescaled so that they average to 1; 
# (iii) the loglikelihood may have multiple local maxima, so ideally a number of different starting values should be tried;
# (iv) when fitting an 'incorrect' model (for example, fitting a model of secondary contact to data generated from a model with decreasing gene flow), 
#      the likelihood maximization function may try smaller and smaller values of the population size parameters, resulting in infinite or NaN values; 
#      if this happens, we suggest increasing the lower bounds on the population size parameters in the optional input vector lb; 
#      trying different starting values may also help.

gim<-function(x1,x2,x3,r1,r2,r3,model,iv=get(paste0("iv.",model)),lb=get(paste0("lb.",model))){
  
  r.average<-mean(c(r1,r2,r3))    # rescale the mutation scalars so that they average to 1
  r1_scaled<-r1/r.average
  r2_scaled<-r2/r.average
  r3_scaled<-r3/r.average

  negll <- get(paste0("negll.",model))
  fit <- nlminb(start=iv,objective=negll,lower=lb,x1=x1,x2=x2,x3=x3,r1=r1_scaled,r2=r2_scaled,r3=r3_scaled,
                control=list(eval.max=1000,iter.max=1000,abs.tol= 1e-20,step.min=0.1,step.max=10,sing.tol=1e-20))
  estimates <- fit$par
  estimates_labelled <- listfunc(model,fit)
  max_lnL <- -fit$objective
  AIC <- 2*(fit$objective+length(estimates))

  printfunc(model,values=estimates)

  cat("\n","maximized loglikelihood =",max_lnL,"\n",
      "AIC score =",AIC,"\n")

  list(fit=fit,estimates=estimates,estimates_labelled=estimates_labelled,max_lnL=max_lnL,AIC=AIC)

}


### Function to compute the Wald standard errors of the ML estimates;
#   this function also computes the covariance matrix.
#
# Required input:
# x1, x2, x3, r1, r2, r3, model: as in the function "gim" above;
# estimates: the vector of (unlabelled) ML estimates obtained with the function "gim" above. 
# 
# The function returns the following:
# se = vector of standard errors of the ML estimates;
# Cov = covariance matrix of the ML estimates.
#
# The function also prints the standard errors of the ML estimates.
# 
# Note: 
# This function can only be used if all the ML estimates are non-zero.
# If one or more parameters are estimated to be zero, then an alternative approach is to reduce the model by fixing at zero the parameter(s) which have 
# been estimated to be zero, and to compute the Wald standard errors for the parameters of the reduced model. As the data will have the same likelihood 
# under both models and the reduced model has fewer parameters, the latter model is preferable.

gim_se<-function(x1,x2,x3,r1,r2,r3,model,estimates){

  H <- gim_hessian(x1,x2,x3,r1,r2,r3,model,params=estimates)
  Cov <- solve(H)
  se <- sqrt(diag(Cov))

  printfunc_se(model,se)

  list(se=se,Cov=Cov)

}


### Function to list and label the ML estimates:
#   (this function is used by the "gim" function above)

listfunc<-function(model,fit){

  if(model=="GIM"){
    par.list<-list(theta0=fit$par[1],theta1=fit$par[2],theta2=fit$par[3],
                   theta1_prime=fit$par[4],theta2_prime=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M1_star=fit$par[8],M2_star=fit$par[9],M1_prime_star=fit$par[10],M2_prime_star=fit$par[11])
  }

  if(model=="IIM"){
    par.list<-list(theta0=fit$par[1],theta1=fit$par[2],theta2=fit$par[3],
                   theta1_prime=fit$par[4],theta2_prime=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M1_star=fit$par[8],M2_star=fit$par[9])
  }
  
  if(model=="SEC"){
    par.list<-list(theta0=fit$par[1],theta1=fit$par[2],theta2=fit$par[3],
                   theta1_prime=fit$par[4],theta2_prime=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M1_prime_star=fit$par[8],M2_prime_star=fit$par[9])
  }
  
  if(model=="ISO"){
    par.list<-list(theta0=fit$par[1],theta1=fit$par[2],theta2=fit$par[3],
                   theta1_prime=fit$par[4],theta2_prime=fit$par[5],T1=fit$par[6],V=fit$par[7])
  }
  
  return(par.list)
  
}


### Function to print the ML estimates:
#   (this function is used by the "gim" function above)

printfunc<-function(model,values){

  if(model=="GIM"){
    cat("\n","GIM model:","\n\n","ML estimate for (theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1*, M_2*, M_1'*, M_2'*)","\n","=",values,"\n")
  }

  if(model=="IIM"){
    cat("\n","IIM model:","\n\n","ML estimate for (theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1*, M_2*)","\n","=",values,"\n")
  }
  
  if(model=="SEC"){
    cat("\n","secondary contact model:","\n\n","ML estimate for (theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1'*, M_2'*)","\n","=",values,"\n")
  }
  
  if(model=="ISO"){
    cat("\n","isolation model:","\n\n","ML estimate for (theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V)","\n","=",values,"\n")
  }
  
}


### Function to compute the Hessian matrix of the negated loglikelihood at any vector of parameter values:
#   (this function is used by the "gim_se" function above)

gim_hessian<-function(x1,x2,x3,r1,r2,r3,model,params){
  
  r.average<-mean(c(r1,r2,r3))    # rescale the mutation scalars so that they average to 1
  r1_scaled<-r1/r.average
  r2_scaled<-r2/r.average
  r3_scaled<-r3/r.average

  negll <- get(paste0("negll.",model))

  H <- hessian(negll,params,x1=x1,x2=x2,x3=x3,r1=r1_scaled,r2=r2_scaled,r3=r3_scaled)

  H  

}


### Function to print the standard errors of the ML estimates:
#   (this function is used by the "gim_se" function above)

printfunc_se<-function(model,se){

  if(model=="GIM"){
    cat("\n","standard errors of theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1*, M_2*, M_1'*, M_2'*:","\n",se,"\n")
  }

  if(model=="IIM"){
    cat("\n","standard errors of theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1*, M_2*:","\n",se,"\n")
  }
  
  if(model=="SEC"){
    cat("\n","standard errors of theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V, M_1'*, M_2'*:","\n",se,"\n")
  }
  
  if(model=="ISO"){
    cat("\n","standard errors of theta_0, theta_1, theta_2, theta_1', theta_2', T_1, V:","\n",se,"\n")
  }
  
}



### COMPUTATION OF THE NEGATED LOGLIKELIHOOD:
#   the functions below compute -ln(L) for the four models considered in the paper.
#   x1, x2, x3, r1, r2, r3 are the input data and relative mutation rates as in the function "gim" above.

# GIM MODEL:
# parameters 1 to 11 correspond to theta0, theta1, theta2, theta1_prime, theta2_prime, T1, V, M1*, M2*, M1_prime*, M2_prime* in the paper.  

negll.GIM<-function(params,x1,x2,x3,r1,r2,r3){

  a.theta<-params[1]
  theta<-params[2]  
  b.theta<-params[3]
  c1.theta<-params[4]
  c2.theta<-params[5]
  T1<-params[6]
  V<-params[7] 
  T0<-T1+V
  phi1<-params[8]/theta
  phi2<-params[9]/b.theta
  phi1prime<-params[10]/c1.theta
  phi2prime<-params[11]/c2.theta

  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c1.theta+phi1prime),0,phi1prime,1/c1.theta)
  Qt[,2]<-c(0,-(1/c2.theta+phi2prime),phi2prime,1/c2.theta)
  Qt[,3]<-c(phi2prime/2,phi1prime/2,-(phi1prime+phi2prime)/2,0)
  Qt[,4]<-c(0,0,0,0)
  H1<-eigen(Qt,symmetric=FALSE)
  G<-t(H1$vectors)
  Ginv<-solve(G)
  nu<-H1$values

  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1/theta+phi1),0,phi1,1/theta)
  Qt2[,2]<-c(0,-(1/b.theta+phi2),phi2,1/b.theta)
  Qt2[,3]<-c(phi2/2,phi1/2,-(phi1+phi2)/2,0)
  Qt2[,4]<-c(0,0,0,0)
  H2<-eigen(Qt2,symmetric=FALSE)
  C<-t(H2$vectors)
  Cinv<-solve(C)
  nu2<-H2$values
  
  pij1<-Ginv%*%diag(exp(nu*T1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*V))%*%C
  
  GG <- -Ginv%*%diag(G[,4])
  CC <- -Cinv%*%diag(C[,4])
  # relevant entries of DD are 1: upper 3x3 submatrix is identity matrix
  
  alpha <- -nu[1:3]
  beta  <- -nu2[1:3]
  gamma <- 1/a.theta

  term1.1<-(alpha/(alpha+R1))*(R1/(alpha+R1))^X1*(1-ppois(X1,T1*(alpha+R1)))  
  term2.1<-(beta/(beta+R1))*(R1/(beta+R1))^X1*exp(beta*T1)*(ppois(X1,T1*(beta+R1))-ppois(X1,T0*(beta+R1)))
  term3.1<-(gamma/(gamma+r1))*(r1/(gamma+r1))^x1*exp(gamma*T0)*ppois(x1,T0*(gamma+r1))

  term1.2<-(alpha/(alpha+R2))*(R2/(alpha+R2))^X2*(1-ppois(X2,T1*(alpha+R2)))
  term2.2<-(beta/(beta+R2))*(R2/(beta+R2))^X2*exp(beta*T1)*(ppois(X2,T1*(beta+R2))-ppois(X2,T0*(beta+R2)))
  term3.2<-(gamma/(gamma+r2))*(r2/(gamma+r2))^x2*exp(gamma*T0)*ppois(x2,T0*(gamma+r2))
  
  term1.3<-(alpha/(alpha+R3))*(R3/(alpha+R3))^X3*(1-ppois(X3,T1*(alpha+R3)))
  term2.3<-(beta/(beta+R3))*(R3/(beta+R3))^X3*exp(beta*T1)*(ppois(X3,T1*(beta+R3))-ppois(X3,T0*(beta+R3)))
  term3.3<-(gamma/(gamma+r3))*(r3/(gamma+r3))^x3*exp(gamma*T0)*ppois(x3,T0*(gamma+r3))

  loglike1<-sum(log(GG[1,1:3]%*%term1.1+
                      pij1[1,1:3]%*%CC[1:3,1:3]%*%term2.1+
                      (1-(pij1%*%pij2)[1,4])*term3.1 ))
  
  loglike2<-sum(log(GG[2,1:3]%*%term1.2+
                      pij1[2,1:3]%*%CC[1:3,1:3]%*%term2.2+
                      (1-(pij1%*%pij2)[2,4])*term3.2 ))
  
  loglike3<-sum(log(GG[3,1:3]%*%term1.3+
                      pij1[3,1:3]%*%CC[1:3,1:3]%*%term2.3+
                      (1-(pij1%*%pij2)[3,4])*term3.3 ))

  -(loglike1+loglike2+loglike3)
  
}

# IIM MODEL (GIM model with M1_prime*=M2_prime*=0): 
# parameters 1 to 9 correspond to theta0, theta1, theta2, theta1_prime, theta2_prime, T1, V, M1*, M2* in the paper.  

negll.IIM<-function(params,x1,x2,x3,r1,r2,r3){
  
  a.theta<-params[1]
  theta<-params[2]  
  b.theta<-params[3]
  c1.theta<-params[4]
  c2.theta<-params[5]
  T1<-params[6]
  V<-params[7] 
  T0<-T1+V
  phi1<-params[8]/theta
  phi2<-params[9]/b.theta
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)

  nu<-c(-1/c1.theta,-1/c2.theta,0,0)
  G<-matrix(c(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(-1,-1,-1,1)),nrow=4,ncol=4)
  Ginv<-matrix(c(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(1,1,1,1)),nrow=4,ncol=4)
  
  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1/theta+phi1),0,phi1,1/theta)
  Qt2[,2]<-c(0,-(1/b.theta+phi2),phi2,1/b.theta)
  Qt2[,3]<-c(phi2/2,phi1/2,-(phi1+phi2)/2,0)
  Qt2[,4]<-c(0,0,0,0)
  H2<-eigen(Qt2,symmetric=FALSE)
  C<-t(H2$vectors)
  Cinv<-solve(C)
  nu2<-H2$values

  pij1<-Ginv%*%diag(exp(nu*T1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*V))%*%C
  
  CC <- -Cinv%*%diag(C[,4])
  # relevant entries of GG and DD are 1 for this model: upper 3x3 submatrix is identity matrix
  
  alpha1 <- 1/c1.theta
  alpha2 <- 1/c2.theta
  beta <- -nu2[1:3]
  gamma <- 1/a.theta

  term1.1<-(alpha1/(alpha1+r1))*(r1/(alpha1+r1))^x1*(1-ppois(x1,T1*(alpha1+r1)))  
  term2.1<-(beta/(beta+R1))*(R1/(beta+R1))^X1*exp(beta*T1)*(ppois(X1,T1*(beta+R1))-ppois(X1,T0*(beta+R1)))
  term3.1<-(gamma/(gamma+r1))*(r1/(gamma+r1))^x1*exp(gamma*T0)*ppois(x1,T0*(gamma+r1))

  term1.2<-(alpha2/(alpha2+r2))*(r2/(alpha2+r2))^x2*(1-ppois(x2,T1*(alpha2+r2)))
  term2.2<-(beta/(beta+R2))*(R2/(beta+R2))^X2*exp(beta*T1)*(ppois(X2,T1*(beta+R2))-ppois(X2,T0*(beta+R2)))
  term3.2<-(gamma/(gamma+r2))*(r2/(gamma+r2))^x2*exp(gamma*T0)*ppois(x2,T0*(gamma+r2))
  
  term2.3<-(beta/(beta+R3))*(R3/(beta+R3))^X3*exp(beta*T1)*(ppois(X3,T1*(beta+R3))-ppois(X3,T0*(beta+R3)))
  term3.3<-(gamma/(gamma+r3))*(r3/(gamma+r3))^x3*exp(gamma*T0)*ppois(x3,T0*(gamma+r3))
    
  loglike1<-sum(log(term1.1+
                      pij1[1,1:3]%*%CC[1:3,1:3]%*%term2.1+
                      (1-(pij1%*%pij2)[1,4])*term3.1 ))
  
  loglike2<-sum(log(term1.2+
                      pij1[2,1:3]%*%CC[1:3,1:3]%*%term2.2+
                      (1-(pij1%*%pij2)[2,4])*term3.2 ))
  
  loglike3<-sum(log(  CC[3,1:3]%*%term2.3+
                      (1-pij2[3,4])*term3.3 ))
    
  -(loglike1+loglike2+loglike3)
  
}

# MODEL OF SECONDARY CONTACT (GIM model with M1*=M2*=0); 
# parameters 1 to 9 correspond to theta0, theta1, theta2, theta1_prime, theta2_prime, T1, V, M1_prime*, M2_prime* in the paper.  

negll.SEC<-function(params,x1,x2,x3,r1,r2,r3){
  
  a.theta<-params[1]
  theta<-params[2]  
  b.theta<-params[3]
  c1.theta<-params[4]
  c2.theta<-params[5]
  T1<-params[6]
  V<-params[7] 
  T0<-T1+V
  phi1prime<-params[8]/c1.theta
  phi2prime<-params[9]/c2.theta
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
  
  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c1.theta+phi1prime),0,phi1prime,1/c1.theta)
  Qt[,2]<-c(0,-(1/c2.theta+phi2prime),phi2prime,1/c2.theta)
  Qt[,3]<-c(phi2prime/2,phi1prime/2,-(phi1prime+phi2prime)/2,0)
  Qt[,4]<-c(0,0,0,0)
  H1<-eigen(Qt,symmetric=FALSE)
  G<-t(H1$vectors)
  Ginv<-solve(G)
  nu<-H1$values
  
  nu2<-c(-1/theta,-1/b.theta,0,0)
  C<-matrix(c(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(-1,-1,-1,1)),nrow=4,ncol=4)
  Cinv<-matrix(c(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(1,1,1,1)),nrow=4,ncol=4)

  pij1<-Ginv%*%diag(exp(nu*T1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*V))%*%C
  
  GG <- -Ginv%*%diag(G[,4])
  # relevant entries of CC and DD are 1 for this model: upper 3x3 submatrix is identity matrix
    
  alpha <- -nu[1:3]
  beta <- -nu2[1:2]
  gamma <- 1/a.theta

  term1.1<-(alpha/(alpha+R1))*(R1/(alpha+R1))^X1*(1-ppois(X1,T1*(alpha+R1)))  
  term2.1<-((beta/(beta+R1[1:2,]))*(R1[1:2,]/(beta+R1[1:2,]))^X1[1:2,]*exp(beta*T1)
            *(ppois(X1[1:2,],T1*(beta+R1[1:2,]))-ppois(X1[1:2,],T0*(beta+R1[1:2,]))))
  term3.1<-(gamma/(gamma+r1))*(r1/(gamma+r1))^x1*exp(gamma*T0)*ppois(x1,T0*(gamma+r1))

  term1.2<-(alpha/(alpha+R2))*(R2/(alpha+R2))^X2*(1-ppois(X2,T1*(alpha+R2)))
  term2.2<-((beta/(beta+R2[1:2,]))*(R2[1:2,]/(beta+R2[1:2,]))^X2[1:2,]*exp(beta*T1)
            *(ppois(X2[1:2,],T1*(beta+R2[1:2,]))-ppois(X2[1:2,],T0*(beta+R2[1:2,]))))
  term3.2<-(gamma/(gamma+r2))*(r2/(gamma+r2))^x2*exp(gamma*T0)*ppois(x2,T0*(gamma+r2))
  
  term1.3<-(alpha/(alpha+R3))*(R3/(alpha+R3))^X3*(1-ppois(X3,T1*(alpha+R3)))
  term2.3<-((beta/(beta+R3[1:2,]))*(R3[1:2,]/(beta+R3[1:2,]))^X3[1:2,]*exp(beta*T1)
            *(ppois(X3[1:2,],T1*(beta+R3[1:2,]))-ppois(X3[1:2,],T0*(beta+R3[1:2,]))))
  term3.3<-(gamma/(gamma+r3))*(r3/(gamma+r3))^x3*exp(gamma*T0)*ppois(x3,T0*(gamma+r3))
  
   
  loglike1<-sum(log(GG[1,1:3]%*%term1.1+
                      pij1[1,1:2]%*%term2.1+
                      (1-(pij1%*%pij2)[1,4])*term3.1 ))
  
  loglike2<-sum(log(GG[2,1:3]%*%term1.2+
                      pij1[2,1:2]%*%term2.2+
                      (1-(pij1%*%pij2)[2,4])*term3.2 ))
  
  loglike3<-sum(log(GG[3,1:3]%*%term1.3+
                      pij1[3,1:2]%*%term2.3+
                      (1-(pij1%*%pij2)[3,4])*term3.3 ))
  
  -(loglike1+loglike2+loglike3)
  
}

# ISOLATION MODEL (GIM model with M1*=M2*=M1_prime*=M2_prime*=0); 
# parameters 1 to 7 correspond to theta0, theta1, theta2, theta1_prime, theta2_prime, T1, V in the paper.  

negll.ISO<-function(params,x1,x2,x3,r1,r2,r3){

  a.theta<-params[1]
  theta<-params[2]  
  b.theta<-params[3]
  c1.theta<-params[4]
  c2.theta<-params[5]
  T1<-params[6]
  V<-params[7] 
  T0<-T1+V

  loglike1<-sum(log((c1.theta*r1)^x1/(1+c1.theta*r1)^(x1+1)*(1-ppois(x1,(1/c1.theta+r1)*T1))+
                      exp(T1*(1/theta-1/c1.theta))*(theta*r1)^x1/(1+theta*r1)^(x1+1)*(ppois(x1,(1/theta+r1)*T1)-ppois(x1,(1/theta+r1)*T0))+
                      exp(T0/a.theta-T1/c1.theta)*(a.theta*r1)^x1/(1+a.theta*r1)^(x1+1)*ppois(x1,(1/a.theta+r1)*T0)*exp(-V/theta)))
  
  loglike2<-sum(log((c2.theta*r2)^x2/(1+c2.theta*r2)^(x2+1)*(1-ppois(x2,(1/c2.theta+r2)*T1))+
                      exp(T1*(1/b.theta-1/c2.theta))*(b.theta*r2)^x2/(1+b.theta*r2)^(x2+1)*(ppois(x2,(1/b.theta+r2)*T1)-ppois(x2,(1/b.theta+r2)*T0))+
                      exp(T0/a.theta-T1/c2.theta)*(a.theta*r2)^x2/(1+a.theta*r2)^(x2+1)*ppois(x2,(1/a.theta+r2)*T0)*exp(-V/b.theta)))
  
  loglike3<-sum(log((a.theta*r3)^x3/(1+a.theta*r3)^(x3+1)*exp(T0/a.theta)*ppois(x3,(1/a.theta+r3)*T0)))
  
  
  -(loglike1+loglike2+loglike3)

}


### Default values used by the "gim" function above:

### Vector of lower bounds for each model:
#   a lower bound > 0 needs to be set for the population size parameters (these may need to be changed to suit the data under consideration); 
#   the lower bounds for the time and migration rate parameters are set to 0.

lb.GIM<-c(rep(0.01,5),rep(0,6))
lb.IIM<-c(rep(0.01,5),rep(0,4))
lb.SEC<-c(rep(0.01,5),rep(0,4))
lb.ISO<-c(rep(0.01,5),rep(0,2))

### Vector of initial values for each model:
#  (these may need to be changed to more appropriate values for the specific data under consideration).

iv.GIM<-c(rep(5,7),rep(0.3,4))
iv.IIM<-c(rep(5,7),rep(0.3,2))
iv.SEC<-c(rep(5,7),rep(0.3,2))
iv.ISO<-rep(5,7)


