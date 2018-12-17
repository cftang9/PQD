# Testing for positive quadrant dependence.

This repository contains R programs for the article, “Testing for positive quadrant dependence.” 
This article has been submitted for publication. 
Prior to using R programs on this repository, please download the main function [EL_PQD_Library.R](https://github.com/cftang9/PQD/blob/master/EL_PQD_Library.R). 

## Three data sets:

### Educational data: 

[EducationData.csv](https://github.com/cftang9/PQD/blob/master/EducationData.csv) 
(Codes: [EducationData.R](https://github.com/cftang9/PQD/blob/master/EducationData.R) )

### Twins Data:  

[TwinsData.csv](https://github.com/cftang9/PQD/blob/master/TwinsData.csv) 
(Codes: [TwinsData.R](https://github.com/cftang9/PQD/blob/master/TwinsData.R) )

### Stocks Data: 

[StockData.csv](https://github.com/cftang9/PQD/blob/master/StockData.csv) 
(Codes: [StockData.R](https://github.com/cftang9/PQD/blob/master/StockData.R) )

## R program for repeating the simulation study: 
[Simulation.R](https://github.com/cftang9/PQD/blob/master/Simulation.R)
(It might take a long time for a single computer to run 10,000 replications. We used a cluster)

### Illustrative examples: 

[IllustrativeExamples.R](https://github.com/cftang9/PQD/blob/master/IllustrativeExamples.R)

### Reproducing the simulation results. 

#### This is a test

## 
For your own data set, please use these commands after naming the data by X and Y:
```R
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
IndvsPQD(X,Y,graph=TRUE)
```
