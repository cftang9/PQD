# Testing for positive quadrant dependence

This repository contains R programs for the article, “Testing for positive quadrant dependence.” 
This article has been submitted for publication. 
Prior to using R programs on this repository, please download the main function [EL_PQD_Library.R](https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R). 

### Illustrative examples: 

For illustration, we generate a random sample from a Clayton copula, with a user-specifiic Kandell's tau, to test for independence versus positive quadrant dependence (PQD). (For large sample sizes, it may take some time to generate the critical values)      
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R"); 
n = 10; tau = 0.2;
X = RV_CopTau(n, tau, Copula="Clayton"); 
IndvsPQD(X=X[,1],Y=X[,2],graph=TRUE); 
```
A scatter plot and a corresponding pseudo-observations plot between `X` and `Y` and will be produce. 
The empirical-likelihood-based test EL and distance-based test KS, CvM, and AD for PQD will be perform with corresponding test statistics, critical values, and p-values. 

You can change the argument `Copula="Calyton"` in the function `RV_CopTau` above into `Copula="Frank"` and `Copula="Gumbel"` and generate random sample from a Frank and Gumbel copulas, respectively.
All the illustration codes and results are included here. 
[IllustrativeExamples.R](https://raw.githubusercontent.com/cftang9/PQD/master/IllustrativeExamples.R)



### Perform the tests for our own data: 
For your own data set, please use these commands after naming the data by X and Y:
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
IndvsPQD(X,Y,graph=TRUE)
```


## Data analysis:

### Educational data: 

[EducationData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.csv)
(Codes: [EducationData.R](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.R) )

### Twins Data:  

[TwinsData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.csv) 
(Codes: [TwinsData.R](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.R) )

### Stocks Data: 

[StockData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.csv) 
(Codes: [StockData.R](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.R) )

## Size and power demonstrations: 

### Clayton, Frank, Gumbel, and Gaussian copulas: 
Clayton, Frank, Gumbel, and Gaussian copulas are classic copulas which can be shown satisfying PQD conditions with non-negative Kandell's tau or non-negative Pearson's rho for Gaussian. 
For each copula, 100,000 Monte Carlo random samples with sample size
[n=100](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D100.R)
were used to estimate the sizes or powers for each test. 
The size and power results with sample sizes 
[n=50](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D50.R)
and 
[n=200](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D200.R)
were provided in the supplementary file. 

### Restricted Students' t-distribution, Farlie–Gumbel–Morgenstern (FGM), and Cuadras–Aug ́e (CA) copulas: 
More copula families were introduced to perform the size and power studies. 
It can be shown that the Students' t-distribution restricted on the first quadrant and centered at the origin with identity shape matrix is PQD. Further, PQD also holds for some FGM and CA copulas. 

Similarly, for each copula, 100,000 Monte Carlo random samples with sample size
[n=100](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D100.R)
were used to estimate the sizes or powers for each test.
The size and power results with sample sizes 
[n=50](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D50.R)
and
[n=200](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D200.R)
were provided in the supplementary file. 

