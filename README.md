# Testing for positive quadrant dependence

This repository contains R programs for the article, “Testing for positive quadrant dependence.” 
This article has been submitted for publication. 
Prior to using R programs on this repository, please download the main R program [EL_PQD_Library.R](https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R) or source it in R using the command: source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R"); 

## Part 1:  Illustration

### 1.1  A simple example: 

Below generates a random sample of size 10 from a Clayton copula, with a user-specified Kendall's tau, to test for independence versus positive quadrant dependence (PQD). 
```R
# Source the main R program
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R"); 
# Set sample size n, and the Kendall's tau
n = 10; tau = 0.2;
# Generate a sample of size n
Sample = RV_CopTau(n, tau, Copula="Clayton"); 
# Name the sample by X and Y
X=Sample[,1];Y=Sample[,2];
# Run the test
IndvsPQD(X,Y,graph=TRUE); 
```
A scatter plot and a plot of the corresponding pseudo-observations between `X` and `Y` and will be produced. 
The empirical-likelihood-based test EL and distance-based test KS, CvM, and AD for PQD will be performed. Results include the corresponding test statistics, critical values at significance level 0.05, and p-values.

The argument `Copula="Calyton"` in the function `RV_CopTau` above can be changed to `Copula="Frank"` and `Copula="Gumbel"` to generate a random sample from the Frank and Gumbel copulas, respectively. The gaussian copula can be considered, too. See these details in [IllustrativeExamples.R](https://raw.githubusercontent.com/cftang9/PQD/master/IllustrativeExamples.R)

For an easy illustration, we set n=10 above. Other sample sizes can be considered as well. However, When the sample size is large, it will take a longer time to run.


### 1.2 For your own data: 
Please use these R commands after naming the data by X and Y:
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
# name your data by X and Y
IndvsPQD(X,Y,graph=TRUE)
```

## Part 2: To repeat the simulation results in Section 3 of the manuscript: 

### 2.1 Clayton, Frank, Gumbel, and Gaussian copulas: 
These are classic copulas. It can be shown that these copulas satisfy the PQD condition if the associated Kandell's tau or Pearson's rho is non-negative.
For each copula, 100,000 Monte Carlo random samples with sample size
[n=100](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D100.R)
were used to estimate the sizes or powers for each test. 
The size and power results with sample sizes 
[n=50](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D50.R)
and 
[n=200](https://raw.githubusercontent.com/cftang9/PQD/master/Clayton%20Frank%20Gumbel%20and%20Gaussian%20n%3D200.R)
were provided in the supplementary file. 

### 2.2 Restricted Students' t-distribution, Farlie–Gumbel–Morgenstern (FGM), and Cuadras–Aug ́e (CA) copulas: 
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


## Part 3: To repeat the real data analysis results in Section 4 of the manuscript:
We applied all tests in this manuscript to three data applications. To repeat the results of our analysis, simply run the R program for each. The data included in the csv file will be automatically read by the corresponding R program.


### 3.1 Twins Data:  

Data: [TwinsData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.csv) 
(R program: [TwinsData.R](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.R) )

### 3.2 Education data: 

Data: [EducationData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.csv)
(R program: [EducationData.R](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.R) )


### 3.3 Stock Data: 

Data: [StockData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.csv) 
(R program: [StockData.R](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.R) )



