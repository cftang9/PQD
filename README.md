# Testing for positive quadrant dependence

This repository contains R programs for the article, “Testing for positive quadrant dependence.” 
This article has been submitted for publication. 
Prior to using R programs on this repository, please download the main function [EL_PQD_Library.R](https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R). 

### Illustrative examples: 

Below generates a random sample of size 10 from a Clayton copula, with a user-specified Kandell's tau, to test for independence versus positive quadrant dependence (PQD). 
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R"); 
n = 10; tau = 0.2;
Sample = RV_CopTau(n, tau, Copula="Clayton"); 
X=Sample[,1];Y=Sample[,2];
IndvsPQD(X,Y,graph=TRUE); 
```
A scatter plot and a plot of the corresponding pseudo-observations between `X` and `Y` and will be produce. 
The empirical-likelihood-based test EL and distance-based test KS, CvM, and AD for PQD will be performed. Results include the corresponding test statistics, critical values at significance level 0.05, and p-values.

You can change the argument `Copula="Calyton"` in the function `RV_CopTau` above to `Copula="Frank"` and `Copula="Gumbel"` and generate a random sample from Frank and Gumbel copulas, respectively. Other sample sizes can be considered as well. However, When the sample size is large, it may take some time to compute the critical values.
All the illustration codes and results are included here. 
[IllustrativeExamples.R](https://raw.githubusercontent.com/cftang9/PQD/master/IllustrativeExamples.R)



### Perform the tests for our own data: 
For your own data set, please use these commands after naming the data by X and Y:
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
# name your data by X and Y
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
These are classic copulas. It can be shown that these copulas satisfy the PQD condition if the associated Kandell's tau or Pearson's rho is non-negative.
For each copula, 100,000 Monte Carlo random samples with sample size
[n=100](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D100.R)
were used to estimate the sizes or powers for each test. 
The size and power results with sample sizes 
[n=50](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D50.R)
and 
[n=200](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D200.R)
were provided in the supplementary file. 

### Restricted Students' t-distribution, Farlie–Gumbel–Morgenstern (FGM), and Cuadras–Aug ́e (CA) copulas: 
More copulas were considered in the size and power studies. 
It can be shown that the Students' t-distribution restricted on the first quadrant and centered at the origin with identity shape matrix satisfies the PQD condition. Further, PQD also holds for some FGM and CA copulas. 

For each of these copulas, 100,000 Monte Carlo random samples with sample size
[n=100](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D100.R)
were used to estimate the sizes or powers for each test.
The size and power results with sample sizes 
[n=50](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D50.R)
and
[n=200](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted%20t%20FGM%20and%20CA%20n%3D200.R)
were provided in the supplementary file. 

