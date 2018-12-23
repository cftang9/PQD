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
    X = rCopula(n, myCop.Clay)
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

Vexler2014 = function(x,y){
  n=length(x); 
  m <- r <- rep(round(0.5*n^0.8), n); i = 1:n; 
  #x = runif(n); y = runif(n); 
  sx = sort(x); sy = sort(y); 
  yD <- 1*((i-m) < 1) + (i-m)*((i-m) >= 1); 
  yU <- n*((i+m) > n) + (i+m)*((i+m) <= n); 
  xt <- x[order(y)]; si <- rank(xt); 
  xD = 1*((si - r) < 1) + (si - r)*((si - r) >= 1)
  xU = n*((si + r) > n) + (si + r)*((si + r) <= n)
  sxxU = sx[xU]; sxxD = sx[xD]; syyU = sy[yU]; syyD = sy[yD]
  
  denom <- xU - xD; 
  denom[denom==0] <- 1/n; 
  # Write a function for calculating the bivariate empirical function based on (Xi, Yi), i=1,...,n #(without n^(-1))
  fn <- function(u, v){
    sum(1*(x< u & y < v) + 0.5*(x == u & y < v) + 0.5*(x < u & y == v) + 0.25*(x == u & y == v))
  }
  # Write a function for the calculation of Equation (2) 
  fnn <- function(s, n){
    delta <- ( fn(sxxU[s] , syyU[s]) - fn(sxxD[s] , syyU[s]) - fn(sxxU[s] , syyD[s]) + fn(sxxD[s] , syyD[s]) + n^0.55)/denom[s]
    if (delta==0) delta <- 1/n^2
    return(delta)
  }
  Vexler2014 <- sum(log((n^0.2)*sapply(i, function(s) {fnn(s, n)})))
  return(Vexler2014)
}

CV_Vexler2014 = function(n,alpha=0.05){
  if(alpha==0.2 ){ri=2}
  if(alpha==0.1 ){ri=3}
  if(alpha==0.05){ri=4}
  if(alpha==0.01){ri=5}
  Table = array(c(
    5,  3.7592,  3.7982,  4.0709,  4.3405, 
    7,  5.7211,  5.9138,  6.0390,  6.3010,
    10,  7.4321,  7.6529,  7.8549,  8.2521,
    15, 10.8383, 11.1268, 11.3766, 11.8930,
    17, 11.6964, 11.9860, 12.2427, 12.7964,
    20, 14.0863, 14.4149, 14.7061, 15.3051,
    23, 15.6590, 15.9941, 16.2962, 16.9399,
    25, 16.6106, 16.9430, 17.2632, 17.9485,
    30, 19.7139, 20.0724, 20.4089, 21.1245,
    35, 22.7744, 23.1606, 23.5258, 24.2849,
    40, 25.7939, 26.1989, 26.5672, 27.3895,
    45, 28.7824, 29.2060, 29.6057, 30.4312,
    50, 32.1148, 32.5648, 32.9714, 33.8318,
    60, 37.9011, 38.3995, 38.8356, 39.7367,
    70, 43.6649, 44.1765, 44.6390, 45.6783,
    80, 49.4025, 49.9159, 50.3733, 51.3542,
    90, 55.3020, 55.8526, 56.3487, 57.4026,
    100, 60.9451, 61.4996, 62.0054, 63.0367
  ), c(5,18))
  ni = c(Table[1,]==n); 
  CV_Vexler2014 = Table[ri,ni]
  return(CV_Vexler2014)
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
