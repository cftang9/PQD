#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
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
# dependence
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
# Function to generate a random sample from Frank, Clayton, or Gumbel copula
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
    GumPara = iTau(gumbelCopula(2),Tau); 
    myCop.Gumbel = archmCopula(family = "gumbel", dim = 2, param = GumPara)
    X = rCopula(n, myCop.Gumbel)
  }
  return(X)
}

#***********************************************************************#
# Function to generate a random sample from the Gaussian copula
# according to Pearson's rho
RV_CopGaussian = function(n=10,Rho=0){
  myCop.Gaussian = normalCopula(Rho)
  X = rCopula(n, myCop.Gaussian)
  return(X)
}

#***********************************************************************#
# Main function to perform the EL, KS, CvM, AD, one-sided Spearman's and 
# Kendall's rank tests
# Arguments:
# X 		   = First part of the paired data. 
# Y 		   = Second part of the paired data. 
# alapha   = Significance level (default 0.05). 
# graph    = TRUE: to produce the scatterplot of the data and the pseudo-observations
#            or FALSE (default): no scatterplot
#
# Value: 	  
# IndvsPQD = A table in data frame provides the test statistic, decision 
#            of rejecting independence (1: reject independence; 0: do not 
#            reject independence), critical value, and p-value for each 
#            following test: EL, KS, CvM, AD, one-sided Spearman and 
#            Kendall's rank tests
# (Large sample size (n=100) may take some time)
IndvsPQD = function(X, Y, alpha=0.05, graph=FALSE){
  n = length(X); 
  if(graph==TRUE)
  {
    par(mar=c(4.5,5,3,0.5))
    par(mfrow=c(1,2))
    plot(X,Y,xlab="X", ylab="Y",xlim=range(X),ylim=range(Y),main="Scatterplot of the data")
    plot(rank(X)/(n+1),rank(Y)/(n+1),xlab="X", ylab="Y",xlim=c(0,1),ylim=c(0,1),main="Scatterplot of pseudo-observations")
  }
  Tn = array(,6); 
  Tn[1] = Rcpp_IP_FF(cbind(X,Y), n); 
  Tn[2:4] = Rcpp_KCA_01(cbind(X,Y), n); 
  Tn[5] = cor(X,Y,method="spearman"); 
  Tn[6] = cor(X,Y,method="kendall"); 
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
  colnames(IndvsPQD) = c("test_statistic", "p-value", "reject_independence", "critical_value")
  rownames(IndvsPQD) = c("EL", "KS", "CvM", "AD", "spearman", "kendall"); 
  print("1: reject independence; 0: do not rejct independence")
  return(IndvsPQD)
}


