### FUNCTIONS TO SIMULATE DATA FROM THE GIM MODEL
### (Costa & Wilkinson-Herbots (2021): 
### "Inference of gene flow in the process of speciation: efficient maximum-likelihood implementation of a generalised isolation-with-migration model".)


### This file contains four versions of a function to simulate data (numbers of nucleotide differences between pairs of sequences) in the GIM model:
### 1) input parameters are the reparameterised ones given by eq.(12) in the paper (see also Figure 2):
###    a) with mutation rate heterogeneity;
###    b) without mutation rate heterogeneity;
### 2) input parameters are the original model parameters (see Figure 1d of the paper):
###    a) with mutation rate heterogeneity;
###    b) without mutation rate heterogeneity.
### These functions can also be used to simulate data from simpler models nested within the GIM model.


### 1a) Function to simulate data (numbers of nucleotide differences between pairs of sequences) in the GIM model with mutation rate heterogeneity 
#      (for each locus, the mutation rate theta is multiplied by a scalar generated from a Gamma distribution):
#
# The input parameters are as in eq.(12) and Figure 2 of the paper;
# nums is a vector with 3 entries; nums[s] is the number of observations to be simulated from locations corresponding to state s 
# (state 1: two sequences from population 1; state 2: two sequences from population 2; state 3: one sequence from each population);
# alpha and lambda are the shape and rate parameters of the Gamma distribution.
#
# The function returns the following: 
# x1, x2, x3: vectors containing the simulated numbers of nucleotide differences between pairs of sequences from states 1, 2 and 3, respectively;
# t1, t2, t3: vectors containing the simulated coalescence times of pairs of sequences from states 1, 2 and 3, respectively;
# r1, r2, r3: vectors containing the simulated mutation scalars for the observations generated from states 1, 2 and 3 respectively;
# n1, n2, n3: numbers of observations simulated from states 1, 2 and 3 respectively.
#
# Note: 
# For simplicity it is preferable for alpha and lambda to be given equal values; 
# in that case the Gamma distribution has mean 1 so that the mutation rate scalars can be interpreted as relative mutation rates.
# If alpha and lambda are not equal then the mutation rate scalars have mean alpha/lambda, so that the mean scaled mutation rate becomes theta*alpha/lambda.

sim_r.GIM<-function(theta0,theta1,theta2,theta1_prime,theta2_prime,T1,V,M1_star,M2_star,M1_prime_star,M2_prime_star,nums,alpha,lambda){ 

  theta = theta1
  a = theta0/theta
  b = theta2/theta
  c1 = theta1_prime/theta
  c2 = theta2_prime/theta
  tau1 = T1/theta
  tau0 = (T1+V)/theta
  M1 = M1_star
  M2 = M2_star/b
  M1c = M1_prime_star/c1
  M2c = M2_prime_star/c2
  
  for (s in 1:3){
    n<-nums[s]
    coal.times<-vector(mode="numeric",length=n)
    for (i in 1:n){ 
      d<-s
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1/c1+M1c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c1)/(1/c1+M1c)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/c2+M2c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c2)/(1/c2+M2c)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1c>0 | M2c>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1c+M2c)/2)
            if(coal.times[i]>tau1){
              coal.times[i]<-tau1
              break
            }
            if(runif(1)<=M2c/(M1c+M2c)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1c==0 & M2c==0){
            coal.times[i]<-tau1
            break
          }
        }
        
      }
      
      if(d==4){next}
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1+M1)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1)/(1+M1)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/b+M2)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1/b)/(1/b+M2)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1>0 | M2>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1+M2)/2)  
            
            if(coal.times[i]>tau0){
              coal.times[i]<-tau0
              break
            }
            if(runif(1)<=M2/(M1+M2)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1==0 & M2==0){
            coal.times[i]<-tau0
            break
          }
        }
        
        
      }
      
      if(d==4){next}else{
        coal.times[i]<-coal.times[i]+rexp(1,1/a)
      }
      
      
    }
    
    r<-rgamma(n,alpha,lambda)
  
    mutations<-rpois(n,theta*r*coal.times)
    if(s==1){
      x1<-mutations
      t1<-coal.times
      n1<-n
      r1<-r
    }
    if(s==2){
      x2<-mutations
      t2<-coal.times
      n2<-n
      r2<-r
    }
    if(s==3){
      x3<-mutations
      t3<-coal.times
      n3<-n
      r3<-r
    }
  }
  
  list(x1=x1,x2=x2,x3=x3,t1=t1,t2=t2,t3=t3,r1=r1,r2=r2,r3=r3,n1=n1,n2=n2,n3=n3)
} 


### 1b) Function to simulate data (numbers of nucleotide differences between pairs of sequences) in the GIM model without mutation rate heterogeneity: 
#
# The input parameters are as in eq.(12) and Figure 2 of the paper;
# nums is a vector with 3 entries; nums[s] is the number of observations to be simulated from locations corresponding to state s 
# (state 1: two sequences from population 1; state 2: two sequences from population 2; state 3: one sequence from each population).
#
# The function returns the following: 
# x1, x2, x3: vectors containing the simulated numbers of nucleotide differences between pairs of sequences from states 1, 2 and 3, respectively;
# t1, t2, t3: vectors containing the simulated coalescence times of pairs of sequences from states 1, 2 and 3, respectively;
# r1, r2, r3: vectors containing the simulated mutation scalars for the observations generated from states 1, 2 and 3 respectively;
# n1, n2, n3: numbers of observations simulated from states 1, 2 and 3 respectively.

sim.GIM<-function(theta0,theta1,theta2,theta1_prime,theta2_prime,T1,V,M1_star,M2_star,M1_prime_star,M2_prime_star,nums){ 

  theta = theta1
  a = theta0/theta
  b = theta2/theta
  c1 = theta1_prime/theta
  c2 = theta2_prime/theta
  tau1 = T1/theta
  tau0 = (T1+V)/theta
  M1 = M1_star
  M2 = M2_star/b
  M1c = M1_prime_star/c1
  M2c = M2_prime_star/c2
  
  for (s in 1:3){
    n<-nums[s]
    coal.times<-vector(mode="numeric",length=n)
    for (i in 1:n){ 
      d<-s
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1/c1+M1c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c1)/(1/c1+M1c)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/c2+M2c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c2)/(1/c2+M2c)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1c>0 | M2c>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1c+M2c)/2)
            if(coal.times[i]>tau1){
              coal.times[i]<-tau1
              break
            }
            if(runif(1)<=M2c/(M1c+M2c)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1c==0 & M2c==0){
            coal.times[i]<-tau1
            break
          }
        }
        
      }
      
      if(d==4){next}
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1+M1)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1)/(1+M1)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/b+M2)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1/b)/(1/b+M2)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1>0 | M2>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1+M2)/2)  
            
            if(coal.times[i]>tau0){
              coal.times[i]<-tau0
              break
            }
            if(runif(1)<=M2/(M1+M2)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1==0 & M2==0){
            coal.times[i]<-tau0
            break
          }
        }
        
        
      }
      
      if(d==4){next}else{
        coal.times[i]<-coal.times[i]+rexp(1,1/a)
      }
      
      
    }
    
    mutations<-rpois(length(coal.times),theta*coal.times)
    if(s==1){
      x1<-mutations
      t1<-coal.times
      n1<-n
    }
    if(s==2){
      x2<-mutations
      t2<-coal.times
      n2<-n
    }
    if(s==3){
      x3<-mutations
      t3<-coal.times
      n3<-n
    }
  }
  
  list(x1=x1,x2=x2,x3=x3,t1=t1,t2=t2,t3=t3,n1=n1,n2=n2,n3=n3)
} 


### 2a) Function to simulate data (numbers of nucleotide differences between pairs of sequences) in the GIM model with mutation rate heterogeneity 
#       (for each locus, the mutation rate theta is multiplied by a scalar generated from a Gamma distribution):
#
# The parameters a, b, c1, c2, tau1, tau0, M1, M2 are as in Figure 1d in the paper; M1c and M2c correspond to M1_prime and M2_prime; 
# theta is the scaled mutation rate;
# nums is a vector with 3 entries; nums[s] is the number of observations to be simulated from locations corresponding to state s 
# (state 1: two sequences from population 1; state 2: two sequences from population 2; state 3: one sequence from each population);
# alpha and lambda are the shape and rate parameters of the Gamma distribution.
#
# The function returns the following: 
# x1, x2, x3: vectors containing the simulated numbers of nucleotide differences between pairs of sequences from states 1, 2 and 3, respectively;
# t1, t2, t3: vectors containing the simulated coalescence times of pairs of sequences from states 1, 2 and 3, respectively;
# r1, r2, r3: vectors containing the simulated mutation scalars for the observations generated from states 1, 2 and 3 respectively;
# n1, n2, n3: numbers of observations simulated from states 1, 2 and 3 respectively.
#
# Note: 
# For simplicity it is preferable for alpha and lambda to be given equal values; 
# in that case the Gamma distribution has mean 1 so that the mutation rate scalars can be interpreted as relative mutation rates.
# If alpha and lambda are not equal then the mutation rate scalars have mean alpha/lambda, so that the mean scaled mutation rate becomes theta*alpha/lambda.

sim_r.GIM.original<-function(a,theta,b,c1,c2,tau1,tau0,M1,M2,M1c,M2c,nums,alpha,lambda){ 
  
  for (s in 1:3){
    n<-nums[s]
    coal.times<-vector(mode="numeric",length=n)
    for (i in 1:n){ 
      d<-s
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1/c1+M1c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c1)/(1/c1+M1c)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/c2+M2c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c2)/(1/c2+M2c)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1c>0 | M2c>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1c+M2c)/2)
            if(coal.times[i]>tau1){
              coal.times[i]<-tau1
              break
            }
            if(runif(1)<=M2c/(M1c+M2c)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1c==0 & M2c==0){
            coal.times[i]<-tau1
            break
          }
        }
        
      }
      
      if(d==4){next}
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1+M1)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1)/(1+M1)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/b+M2)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1/b)/(1/b+M2)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1>0 | M2>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1+M2)/2)  
            
            if(coal.times[i]>tau0){
              coal.times[i]<-tau0
              break
            }
            if(runif(1)<=M2/(M1+M2)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1==0 & M2==0){
            coal.times[i]<-tau0
            break
          }
        }
        
        
      }
      
      if(d==4){next}else{
        coal.times[i]<-coal.times[i]+rexp(1,1/a)
      }
      
      
    }
    
    r<-rgamma(n,alpha,lambda)
  
    mutations<-rpois(n,theta*r*coal.times)
    if(s==1){
      x1<-mutations
      t1<-coal.times
      n1<-n
      r1<-r
    }
    if(s==2){
      x2<-mutations
      t2<-coal.times
      n2<-n
      r2<-r
    }
    if(s==3){
      x3<-mutations
      t3<-coal.times
      n3<-n
      r3<-r
    }
  }
  
  list(x1=x1,x2=x2,x3=x3,t1=t1,t2=t2,t3=t3,r1=r1,r2=r2,r3=r3,n1=n1,n2=n2,n3=n3)
} 


### 2b) Function to simulate data (numbers of nucleotide differences between pairs of sequences) in the GIM model without mutation rate heterogeneity:
#
# The parameters a, b, c1, c2, tau1, tau0, M1, M2 are as in Figure 1d in the paper; M1c and M2c correspond to M1_prime and M2_prime; 
# theta is the scaled mutation rate;
# nums is a vector with 3 entries; nums[s] is the number of observations to be simulated from locations corresponding to state s 
# (state 1: two sequences from population 1; state 2: two sequences from population 2; state 3: one sequence from each population).
#
# The function returns the following: 
# x1, x2, x3: vectors containing the simulated numbers of nucleotide differences between pairs of sequences from states 1, 2 and 3, respectively;
# t1, t2, t3: vectors containing the simulated coalescence times of pairs of sequences from states 1, 2 and 3, respectively;
# n1, n2, n3: numbers of observations simulated from states 1, 2 and 3 respectively.

sim.GIM.original<-function(a,theta,b,c1,c2,tau1,tau0,M1,M2,M1c,M2c,nums){ 
  
  for (s in 1:3){
    n<-nums[s]
    coal.times<-vector(mode="numeric",length=n)
    for (i in 1:n){ 
      d<-s
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1/c1+M1c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c1)/(1/c1+M1c)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/c2+M2c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c2)/(1/c2+M2c)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1c>0 | M2c>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1c+M2c)/2)
            if(coal.times[i]>tau1){
              coal.times[i]<-tau1
              break
            }
            if(runif(1)<=M2c/(M1c+M2c)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1c==0 & M2c==0){
            coal.times[i]<-tau1
            break
          }
        }
        
      }
      
      if(d==4){next}
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1+M1)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1)/(1+M1)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/b+M2)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1/b)/(1/b+M2)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1>0 | M2>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1+M2)/2)  
            
            if(coal.times[i]>tau0){
              coal.times[i]<-tau0
              break
            }
            if(runif(1)<=M2/(M1+M2)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1==0 & M2==0){
            coal.times[i]<-tau0
            break
          }
        }
        
        
      }
      
      if(d==4){next}else{
        coal.times[i]<-coal.times[i]+rexp(1,1/a)
      }
      
      
    }
    
    mutations<-rpois(length(coal.times),theta*coal.times)
    if(s==1){
      x1<-mutations
      t1<-coal.times
      n1<-n
    }
    if(s==2){
      x2<-mutations
      t2<-coal.times
      n2<-n
    }
    if(s==3){
      x3<-mutations
      t3<-coal.times
      n3<-n
    }
  }
  
  list(x1=x1,x2=x2,x3=x3,t1=t1,t2=t2,t3=t3,n1=n1,n2=n2,n3=n3)
} 




