library(Rcpp)
library(copula)

#***********************************************************************#
# Function to produce EL-based test statistic for independence versus 
# positive quadrant dependence. 
cppFunction(" double Rcpp_IP_FF(NumericMatrix X, int n){
            NumericVector Pn11(n*n);
            NumericVector Pn12(n*n);
            NumericVector Pn21(n*n);
            NumericVector Pn22(n*n);
            NumericVector Fn(n*n);
            NumericVector Gn(n*n);
            double Tn=0;
            int ii; int ij; 
            for(int i=0; i<n*n; i++){
            ii = i/n; ij = i%n; 
            for(int j=0; j<n; j++){
            if( X(j,0)<= X(ii,0) & X(j,1)<= X(ij,1)){
            Pn11(i) = Pn11(i)+1; 
            }
            if( X(j,0)<= X(ii,0) & X(j,1)>  X(ij,1)){
            Pn12(i) = Pn12(i)+1; 
            }  
            if( X(j,0)>  X(ii,0) & X(j,1)<= X(ij,1)){
            Pn21(i) = Pn21(i)+1; 
            }  
            if( X(j,0)>  X(ii,0) & X(j,1)>  X(ij,1)){
            Pn22(i) = Pn22(i)+1; 
            }  
            }
            Pn11(i) = Pn11(i)/n; Pn12(i) = Pn12(i)/n; 
            Pn21(i) = Pn21(i)/n; Pn22(i) = Pn22(i)/n;
            Fn(i) = Pn11(i)+Pn12(i); Gn(i) = Pn11(i)+Pn21(i);
            if(Pn11(i)>Fn(i)*Gn(i)){
            if(Pn11(i)*Fn(i)*Gn(i)>0){
            Tn = Tn + n*Pn11(i)*log(Fn(i)*Gn(i)/Pn11(i));
            }
            if(Pn12(i)*Fn(i)*(1-Gn(i))>0){
            Tn = Tn + n*Pn12(i)*log(Fn(i)*(1-Gn(i))/Pn12(i));
            }
            if(Pn21(i)*(1-Fn(i))*Gn(i)>0){
            Tn = Tn + n*Pn21(i)*log((1-Fn(i))*Gn(i)/Pn21(i));
            }
            if(Pn22(i)*(1-Fn(i))*(1-Gn(i))>0){
            Tn = Tn + n*Pn22(i)*log((1-Fn(i))*(1-Gn(i))/Pn22(i));
            }
            }
            }
            return -2*Tn/(n*n);
            }")

#***********************************************************************#
# Function to produce EL-based test statistic of omnibus independence 
# test proposed in Einmahl and McKeague (2003)
cppFunction(" double Rcpp_Ind_FF(NumericMatrix X, int n){
            NumericVector Pn11(n*n);
            NumericVector Pn12(n*n);
            NumericVector Pn21(n*n);
            NumericVector Pn22(n*n);
            NumericVector Fn(n*n);
            NumericVector Gn(n*n);
            double Tn=0;
            int ii; int ij; 
            for(int i=0; i<n*n; i++){
            ii = i/n; ij = i%n; 
            for(int j=0; j<n; j++){
            if( X(j,0)<= X(ii,0) & X(j,1)<= X(ij,1)){
            Pn11(i) = Pn11(i)+1; 
            }
            if( X(j,0)<= X(ii,0) & X(j,1)>  X(ij,1)){
            Pn12(i) = Pn12(i)+1; 
            }  
            if( X(j,0)>  X(ii,0) & X(j,1)<= X(ij,1)){
            Pn21(i) = Pn21(i)+1; 
            }  
            if( X(j,0)>  X(ii,0) & X(j,1)>  X(ij,1)){
            Pn22(i) = Pn22(i)+1; 
            }  
            }
            Pn11(i) = Pn11(i)/n; Pn12(i) = Pn12(i)/n; 
            Pn21(i) = Pn21(i)/n; Pn22(i) = Pn22(i)/n;
            Fn(i) = Pn11(i)+Pn12(i); Gn(i) = Pn11(i)+Pn21(i);
            if(Pn11(i)*Fn(i)*Gn(i)>0){
            Tn = Tn + n*Pn11(i)*log(Fn(i)*Gn(i)/Pn11(i));
            }
            if(Pn12(i)*Fn(i)*(1-Gn(i))>0){
            Tn = Tn + n*Pn12(i)*log(Fn(i)*(1-Gn(i))/Pn12(i));
            }
            if(Pn21(i)*(1-Fn(i))*Gn(i)>0){
            Tn = Tn + n*Pn21(i)*log((1-Fn(i))*Gn(i)/Pn21(i));
            }
            if(Pn22(i)*(1-Fn(i))*(1-Gn(i))>0){
            Tn = Tn + n*Pn22(i)*log((1-Fn(i))*(1-Gn(i))/Pn22(i));
            }
            }
            return -2*Tn/(n*n);
            }")

#***********************************************************************#
# Function to produce the Kolmogorov–Smirnov, Cramer-von-Mises, and 
# Anderson-Darling functional distance between empirical copula and 
# independence copula for testing independence versus positive quadrant 
# dependence. 
cppFunction(" NumericVector Rcpp_KCA_01(NumericMatrix X, int n){
            NumericVector Fn(n);
            NumericVector Gn(n);
            for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
            if( X(j,0)<= X(i,0)){
            Fn(i) = Fn(i)+(double)1/(n+1); 
            }
            if( X(j,1)<= X(i,1)){
            Gn(i) = Gn(i)+(double)1/(n+1); 
            }
            }
            }
            NumericVector Cn(n);
            for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
            if( Fn(j) <= Fn(i) & Gn(j) <= Gn(i)){
            Cn(i) = Cn(i) + (double)1/n; 
            }
            }
            }
            double KS = 0; double CvM = 0; double AD = 0; 
            for(int i=0; i<n; i++){
            if(Cn(i)>Fn(i)*Gn(i)){
            KS = std::max( pow(n,0.5)*(Cn(i)-Fn(i)*Gn(i)), KS);
            CvM = pow(Cn(i)-Fn(i)*Gn(i),2) + CvM;
            AD = pow(Cn(i)-Fn(i)*Gn(i),2)/(Fn(i)*Gn(i)*(1-Fn(i))*(1-Gn(i))) + AD;
            }
            }
            NumericVector TTS(3); 
            TTS(0) = KS; TTS(1) = CvM; TTS(2) = AD;
            return TTS;
            }")

#***********************************************************************#
# Function to generating a random vector sample from Frank and Clayton 
# according to Kandell's tau 
RV_CopTau = function(n=10,Tau=0,Copula="Frank"){
  if(Copula == "Frank"){
    FrankPara = iTau(frankCopula(2),Tau); 
    myCop.Frank = archmCopula(family = "frank", dim = 2, param = FrankPara)
    X = rCopula(n, myCop.Frank)
  }
  if(Copula == "Clayton"){
    ClayPara = iTau(claytonCopula(2),Tau); 
    myCop.Clay = archmCopula(family = "clayton", dim = 2, param = ClayPara)
    X = rCopula(n, myCop.Clay)
  }
  if(Copula == "Gumbel"){
    GumbelPara = iTau(gumbelCopula(2),Tau); 
    myCop.Gumbel = archmCopula(family = "gumbel", dim = 2, param = GumbelPara)
    X = rCopula(n, myCop.Gumbel)
  }
  return(X)
}

#***********************************************************************#
# Function to generating a random vector sample from Gaussian
# according to Pearson's rho
RV_CopGaussian = function(n=10,Rho=0){
  myCop.Gaussian = normalCopula(Rho)
  X = rCopula(n, myCop.Gaussian)
  return(X)
}

#***********************************************************************#
# Main function to perform the EL-based test
# Arguments:
# X 		   = N by 2 data matrix. 
# n        = Number of observations. 
# alapha   = Significance level (defult 0.05). 
#
# Value: 	  
# EL 		   = A list of three variables in data frame. They are the 
#            EL-based test statistic, p-value, and critical value. 
# (Large sample size (n=100) may take some time)
IndvsPQD.EL = function(X, n, alpha=0.05, GenCV=FALSE){
  Tn = Rcpp_IP_FF(X, n); 
  if(GenCV==TRUE){
    B = 10000; Tn0 = array(,B); 
    for(b in 1:B){
      X0 = array(runif(2*n),c(n,2)); 
      Tn0[b] = Rcpp_IP_FF(X0, n); 
    }
    cv_EL = quantile(Tn0,1-alpha); p_EL = mean(Tn0>Tn); 
    EL = data.frame(Tn,p_EL,cv_EL); 
    colnames(EL) = c("test statistic", "p-value", "critical value")
    rownames(EL) = c(""); 
    return(list(EL=EL))
  }
  if(GenCV==FALSE){
    EL = data.frame(Tn); 
    colnames(EL) = c("test statistic")
    rownames(EL) = c(""); 
    return(list(EL=EL))
  }
}

#***********************************************************************#
# Main function to perform the distance-based test
# Arguments:
# X 		   = N by 2 data matrix. 
# n        = Number of observations. 
# alpha    = Significance level
# GenCV    = Generating critical values and p-values
#
# Value: 	  
# KS 		   = A list of variables. If GenCV is TRUE, then a list of three
#            variables: distance-based Kolmogorov–Smirnov test statistic, 
#            p-value, and generated critical value are included. 
#            Otherwise, only the test statistic will be reported. 
# CvM 	   = A list of variables. If GenCV is TRUE, then a list of three
#            variables: distance-based Cramer-von-Mises test statistic, 
#            p-value, and generated critical value are included. 
#            Otherwise, only the test statistic will be reported. 
# AD 		   = A list of variables. If GenCV is TRUE, then a list of three
#            variables: distance-based Anderson-Darling test statistic, 
#            p-value, and generated critical value are included. 
#            Otherwise, only the test statistic will be reported. 
IndvsPQD.DS = function(X, n, alpha=0.05, GenCV=FALSE){
  Tn = Rcpp_KCA_01(X, n); 
  if(GenCV==TRUE){
    B = 10000; Tn0 = array(,c(B,3)); 
    for(b in 1:B){
      X0 = array(runif(2*n),c(n,2)); 
      Tn0[b,] = Rcpp_KCA_01(X0, n); 
    }
    cv_EL = apply(Tn0,2,quantile,prob=1-alpha); p_EL = array(,3); 
    for(j in 1:3){p_EL[j] = mean(Tn0[,j]>Tn[j])}
    KS = data.frame(Tn[1], p_EL[1], cv_EL[1]); 
    colnames(KS)=c("test statistic", "p-value", "critical value")
    rownames(KS) = c("")
    CvM = data.frame(Tn[2], p_EL[2], cv_EL[2]); 
    colnames(CvM)=c("test statistic", "p-value", "critical value")
    rownames(CvM) = c("")
    AD = data.frame(Tn[3], p_EL[3], cv_EL[3]); 
    colnames(AD)=c("test statistic", "p-value", "critical value")
    rownames(AD) = c("")
    return(list(KS = KS, CvM = CvM, AD = AD))
  }
  if(GenCV==FALSE){
    KS = data.frame(Tn[1]); 
    colnames(KS)=c("test statistic")
    rownames(KS) = c("")
    CvM = data.frame(Tn[2]); 
    colnames(CvM)=c("test statistic")
    rownames(CvM) = c("")
    AD = data.frame(Tn[3]); 
    colnames(AD)=c("test statistic")
    rownames(AD) = c("")
    return(list(KS = KS, CvM = CvM, AD = AD))
  }
}

#***********************************************************************#
# Main function to perform the distance-based test
# Arguments:
# X 		   = N by 2 data matrix. 
# n        = Number of observations. 
# alpha    = Significance level
# GenCV    = Generating critical values and p-values
#
# Value: 	  
# spearman = A list of variables. If GenCV is TRUE, then a list of three
#            variables: Spearman's coefficient rho, p-value, and 
#            generated critical value are included. 
#            Otherwise, only the Spearman's coefficient rho will be 
#            reported. 
# kendall  = A list of variables. If GenCV is TRUE, then a list of three
#            variables: Kendall's coefficient rho, p-value, and 
#            generated critical value are included. 
#            Otherwise, only the Kendall's coefficient rho will be 
#            reported. 
IndvsPQD.CO = function(X, n, alpha=0.05, GenCV=FALSE){
  Tn = array(,2)
  Tn[1] = cor(X[,1],X[,2],method="spearman");
  Tn[2] = cor(X[,1],X[,2],method="kendall");
  if(GenCV==TRUE){
    B = 10000; Tn0 = array(,c(B,2)); 
    for(b in 1:B){
      X0 = array(runif(2*n),c(n,2)); 
      Tn0[b,1] = cor(X0[,1],X0[,2],method="spearman");
      Tn0[b,2] = cor(X0[,1],X0[,2],method="kendall");
    }
    cv_EL = apply(Tn0,2,quantile,prob=1-alpha); p_EL = array(,2); 
    for(j in 1:2){p_EL[j] = mean(Tn0[,j]>Tn[j])}
    spearman = data.frame(Tn[1], p_EL[1], cv_EL[1]); 
    colnames(spearman)=c("test statistic", "p-value", "critical value")
    rownames(spearman) = c("")
    kendall = data.frame(Tn[2], p_EL[2], cv_EL[2]); 
    colnames(kendall)=c("test statistic", "p-value", "critical value")
    rownames(kendall) = c("")
    return(list(spearman = spearman, kendall = kendall))
  }
  if(GenCV==FALSE){
    spearman = data.frame(Tn[1]); 
    colnames(spearman)=c("test statistic")
    rownames(spearman) = c("")
    kendall = data.frame(Tn[2]); 
    colnames(kendall)=c("test statistic")
    rownames(kendall) = c("")
    return(list(spearman = spearman, kendall = kendall))
  }
}

IndvsPQD = function(X, Y, CV=NULL, alpha=0.05, graph=FALSE){
  n = length(X); 
  Tn = array(,6); 
  Tn[1] = Rcpp_IP_FF(cbind(X,Y), n); 
  Tn[2:4] = Rcpp_KCA_01(cbind(X,Y), n); 
  Tn[5] = cor(X,Y,method="spearman"); 
  Tn[6] = cor(X,Y,method="kendall"); 
  if(graph==TRUE){
    par(mar=c(4.5,5,3,0.5))
    par(mfrow=c(1,2))
    plot(X, Y, xlab="", ylab="", main="Scatterplot of the data") 
    plot(rank(X)/(n+1),rank(Y)/(n+1), xlab="", ylab="", main="Scatterplot of pseudo-observations"); 
    par(mfrow=c(1,1))
  }
  if(is.null(CV)){
    B = 10000; Tn0 = array(,c(B,6)); 
    for(b in 1:B){
      X0 = array(runif(2*n),c(n,2)); 
      Tn0[b,1] = Rcpp_IP_FF(X0, n); 
      Tn0[b,2:4] = Rcpp_KCA_01(X0, n); 
      Tn0[b,5] = cor(X0[,1],X0[,2],method="spearman"); 
      Tn0[b,6] = cor(X0[,1],X0[,2],method="kendall");
    }
    
    CV = array(,6); p_values = array(,6); Reject_Indep = array(,6)
    CV = apply(Tn0,2,quantile,1-alpha); 
    for(j in 1:6){
      p_values[j] = mean(Tn0[,j]>Tn[j])
      Reject_Indep[j] = (Tn[j]>CV[j])
    }
    IndvsPQD = data.frame(cbind(Tn,p_values,Reject_Indep,CV))
    colnames(IndvsPQD) = c("test statistic", "p-value", "reject independence", "critical value")
    rownames(IndvsPQD) = c("EL", "KS", "CvM", "AD", "spearman", "kendall"); 
    return(IndvsPQD)
  }
  else{
    Reject_Indep = array(,6)
    for(j in 1:6){Reject_Indep[j] = as.logic(as.character(Tn[j]>CV[j]))}
    IndvsPQD = data.frame(cbind(Tn,Reject_Indep,CV))
    colnames(IndvsPQD) = c("test statistic", "reject independence", "critical value")
    rownames(IndvsPQD) = c("EL", "KS", "CvM", "AD", "spearman", "kendall"); 
    return(IndvsPQD)
  }
}

FlipQ = function(X){
  nX = length(X[,1]); 
  for(i in 1:nX){
  if(X[i,1]>1/2 & X[i,2]>1/2){
    X[i,1]=X[i,1]; X[i,2]=X[i,2]; 
  }
  if(X[i,1]<=1/2 & X[i,2]>1/2){
    X[i,1]=1-X[i,1]; X[i,2]=X[i,2]; 
  }
  if(X[i,1]>1/2 & X[i,2]<=1/2){
    X[i,1]=X[i,1]; X[i,2]=1-X[i,2]; 
  }
  if(X[i,1]<=1/2 & X[i,2]<=1/2){
    X[i,1]=1-X[i,1]; X[i,2]=1-X[i,2]; 
  }
  }
  return(X)
}
