rm(list=ls(all=TRUE))
library(xtable)
#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
library(Rcpp)
library(copula)
library(FDGcopulas)
library(GFGM.copula)

source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")

#########################################################################
#  Generating critical values for EL, KS, CvM, AD, Spearman, and Kendall
#  tests for independence versus PQD from uniform (0,1)^2 with sample 
#  size n and 10000 Monte Carlo replications. 
#  (It may take some time to generate critical for large sample sizes, 
#   for example, 100.)
n = 100; #sample size; 
alpha = 0.05; #significance level 

#CV = c(1.425291, 0.6719852, 0.08987011, 4.724722, 0.2342857, 0.1591837) #n=50 set.seed(50)
CV = c(1.402599, 0.6740648, 0.07325512, 3.728709, 0.1646115, 0.1119192) #n=100 set.seed(100)
#CV = c(1.393776, 0.6801778, 0.06367364, 2.948650, 0.1178105, 0.07880402) #n=200 set.seed(200)

#########################################################################
#  Perform simulation using conditional 2D Student's t distribution on the first quadrant
#  with sample size n and B Monte Carlo replications. 
B = 10000#number of Monte Carlo replications
dt = c(1,2,3,4,5,10,50,100); ndt = length(dt); 
PowerT = array(0,c(6,ndt)); 
set.seed(1600)
for(i in 1:ndt){
  myCop.t = tCopula(0,dim=2,df=dt[i]); 
  for(b in 1:B){
    X = rCopula(n, myCop.t); X = FlipQ(X); 
    PowerT[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + PowerT[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    PowerT[2,i] = c(temp[1]>CV[2])/B + PowerT[2,i]; 
    PowerT[3,i] = c(temp[2]>CV[3])/B + PowerT[3,i]; 
    PowerT[4,i] = c(temp[3]>CV[4])/B + PowerT[4,i]; 
    PowerT[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + PowerT[5,i]
    PowerT[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + PowerT[6,i]
  }
}
colnames(PowerT) = dt; rownames(PowerT) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")

#########################################################################
#  Perform simulation using conditional 2D FGM copula on the first quadrant
#  with sample size n and B Monte Carlo replications. 
B = 10000#number of Monte Carlo replications
dt = seq(0,1,by=0.2); ndt = length(dt); 
PowerFGM = array(0,c(6,ndt)); 
set.seed(1700)
for(i in 1:ndt){
  fgm.cop <- fgmCopula(param=dt[i],dim=2)
  for(b in 1:B){
    X <- rCopula(n, fgm.cop)
    PowerFGM[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + PowerFGM[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    PowerFGM[2,i] = c(temp[1]>CV[2])/B + PowerFGM[2,i]; 
    PowerFGM[3,i] = c(temp[2]>CV[3])/B + PowerFGM[3,i]; 
    PowerFGM[4,i] = c(temp[3]>CV[4])/B + PowerFGM[4,i]; 
    PowerFGM[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + PowerFGM[5,i]
    PowerFGM[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + PowerFGM[6,i]
  
  }
}
colnames(PowerFGM) = dt; rownames(PowerFGM) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")

#########################################################################
#  Perform simulation using Cuadras-Aug'e copula distribution on the first quadrant
#  with sample size n and B Monte Carlo replications. 
B = 10000#number of Monte Carlo replications
dt = seq(0,1,by=0.2); ndt = length(dt); 
PowerCA = array(0,c(6,ndt)); 
set.seed(1800)
for(i in 1:ndt){
  theta = c(dt[i], dt[i])
  myFDGcopula <- FDGcopula("cuadrasauge", theta)
  for(b in 1:B){
    X = rFDG(n, myFDGcopula)
    PowerCA[1,i] = c(Rcpp_IP_FF(X,n)>CV[1])/B + PowerCA[1,i]; 
    temp = Rcpp_KCA_01(X,n);  
    PowerCA[2,i] = c(temp[1]>CV[2])/B + PowerCA[2,i]; 
    PowerCA[3,i] = c(temp[2]>CV[3])/B + PowerCA[3,i]; 
    PowerCA[4,i] = c(temp[3]>CV[4])/B + PowerCA[4,i]; 
    PowerCA[5,i] = c(cor(X,method="spearman")[2]>CV[5])/B + PowerCA[5,i]
    PowerCA[6,i] = c(cor(X,method="kendall")[2]>CV[6])/B + PowerCA[6,i]
    #PowerCA[7,i] = c(Vexler2014(X[,1], X[,2]) > CV_Vexler2014(n))/B + PowerCA[7,i]
  }
}
colnames(PowerCA) = dt; rownames(PowerCA) = c("EL","KS","CvM","AD","spearman_rho","kendall_tau")

print("Restricted Students' t"); print(PowerT,digits=3)
print("FGM"); print(PowerFGM,digits=3)
print("CA"); print(PowerCA,digits=3)
