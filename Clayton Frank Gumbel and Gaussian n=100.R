#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
library(Rcpp)
library(copula)
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")

#########################################################################
#  Generating critical values for EL, KS, CvM, AD, Spearman, and Kendall
#  tests for independence versus PQD from uniform (0,1)^2 with sample 
#  size n and 10000 Monte Carlo replications. 
#  (It may take some time to generate critical for large sample sizes, 
#   for example, 100.)
set.seed(100)
n = 100; #sample size; 
alpha = 0.05; #significance level 
CV = array(0,6); TestStats = array(,c(10000,6)); 
for(b in 1:10000){
  X0 = array(runif(2*n),c(n,2)); 
  TestStats[b,1] = Rcpp_IP_FF(X0,n); 
  TestStats[b,2:4] = Rcpp_KCA_01(X0,n);  
  TestStats[b,5] = cor(X0,method="spearman")[2]
  TestStats[b,6] = cor(X0,method="kendall")[2]
}
CV = apply(TestStats,2,quantile,prob=0.95)
CV = data.frame(t(CV))
colnames(CV) = c("EL","KS","CvM","AD","spearman","kendall")

#########################################################################
#  Perform simulation using Clayton copula 
#  with sample size n and B Monte Carlo replications. 
B = 10; #number of Monte Carlo replications
Tau = seq(0,0.4,by=0.1); nTau = length(Tau); 
Power = array(0,c(6,nTau)); 
#set.seed(200) # for n not equals to 200
set.seed(100) # for n=100
for(i in 1:nTau){
  for(b in 1:B){
    if(Tau[i]==0){X = array(runif(2*n),c(n,2)); }
    if(Tau[i]> 0){X = RV_CopTau(n, Tau[i], Copula="Clayton"); }
    Power[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + Power[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    Power[2,i] = c(temp[1]>CV[2])/B + Power[2,i]; 
    Power[3,i] = c(temp[2]>CV[3])/B + Power[3,i]; 
    Power[4,i] = c(temp[3]>CV[4])/B + Power[4,i]; 
    Power[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + Power[5,i]
    Power[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + Power[6,i]
  }
}
colnames(Power) = Tau; rownames(Power) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")
Power_Calyton = Power

#########################################################################
#  Perform simulation using Frank copula 
#  with sample size n and B Monte Carlo replications. 
B = 10; #number of Monte Carlo replications
Tau = seq(0,0.4,by=0.1); nTau = length(Tau); 
Power = array(0,c(6,nTau)); 
set.seed(300)
for(i in 1:nTau){
  for(b in 1:B){
    if(Tau[i]==0){X = array(runif(2*n),c(n,2)); }
    if(Tau[i]> 0){X = RV_CopTau(n, Tau[i], Copula="Frank"); }
    Power[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + Power[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    Power[2,i] = c(temp[1]>CV[2])/B + Power[2,i]; 
    Power[3,i] = c(temp[2]>CV[3])/B + Power[3,i]; 
    Power[4,i] = c(temp[3]>CV[4])/B + Power[4,i]; 
    Power[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + Power[5,i]
    Power[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + Power[6,i]
  }
}
colnames(Power) = Tau; rownames(Power) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")
Power_Frank = Power

#########################################################################
#  Perform simulation using Gumbel copula 
#  with sample size n and B Monte Carlo replications. 
B = 10; #number of Monte Carlo replications
Tau = seq(0,0.4,by=0.1); nTau = length(Tau); 
Power = array(0,c(6,nTau)); 
set.seed(400)
for(i in 1:nTau){
  for(b in 1:B){
    if(Tau[i]==0){X = array(runif(2*n),c(n,2)); }
    if(Tau[i]> 0){X = RV_CopTau(n, Tau[i], Copula="Gumbel"); }
    Power[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + Power[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    Power[2,i] = c(temp[1]>CV[2])/B + Power[2,i]; 
    Power[3,i] = c(temp[2]>CV[3])/B + Power[3,i]; 
    Power[4,i] = c(temp[3]>CV[4])/B + Power[4,i]; 
    Power[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + Power[5,i]
    Power[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + Power[6,i]
  }
}
colnames(Power) = Tau; rownames(Power) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")
Power_Gumbel = Power

#########################################################################
#  Perform simulation using Gaussian copula 
#  with sample size n and B Monte Carlo replications. 
B = 10; #number of Monte Carlo replications
Rho = seq(0,0.4,by=0.1); nRho = length(Rho); 
Power = array(0,c(6,nRho)); 
set.seed(500)
for(i in 1:nRho){
  for(b in 1:B){
    if(Rho[i]==0){X = array(runif(2*n),c(n,2)); }
    if(Rho[i]> 0){X = RV_CopGaussian(n,Rho[i]) }
    Power[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + Power[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    Power[2,i] = c(temp[1]>CV[2])/B + Power[2,i]; 
    Power[3,i] = c(temp[2]>CV[3])/B + Power[3,i]; 
    Power[4,i] = c(temp[3]>CV[4])/B + Power[4,i]; 
    Power[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + Power[5,i]
    Power[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + Power[6,i]
  }
}
colnames(Power) = Tau; rownames(Power) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")
Power_Gaussian = Power




print(Power_Calyton,digit=3); 
print(Power_Frank,digit=3);
print(Power_Gumbel,digit=3);
print(Power_Gaussian,digit=3);


