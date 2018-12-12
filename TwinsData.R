#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
source("http://people.stat.sc.edu/wang528/PQD/EL_PQD_Library.R")
###########################################################################
# This twins data is from Ashenfelter and Krueger (1994). 
raw_data = read.csv("http://people.stat.sc.edu/wang528/PQD/TwinsData.csv")
TwinsData = data.frame(cbind(raw_data$LHRWAGEL, raw_data$LHRWAGEH))
TwinsData = na.omit(TwinsData); 
colnames(TwinsData) = c("LogWage1", "LogWage2")
n = length(TwinsData$LogWage1); #sample size n=149

###########################################################################
# Here we provide the Scatterplot of the data and pseudo-observations 
# to roughly visualize the dependence structure between log wages
par(mar=c(4.5,5,3,0.5))
par(mfrow=c(1,2))
plot(TwinsData$LogWage1, TwinsData$LogWage2, xlab="Twin 1", ylab="Twin 2",xlim=range(TwinsData$LogWage1),ylim=range(TwinsData$LogWage2),main="Scatterplot of the data")
plot(rank(TwinsData$LogWage1)/(n+1), rank(TwinsData$LogWage2)/(n+1), xlab="Twin 1", ylab="Twin 2",xlim=c(0,1),ylim=c(0,1),main="Scatterplot of pseudo-observations")

# Perform all considered tests: 
# It takes <= 10 mins
set.seed(100)
IndvsPQD(TwinsData$LogWage1,TwinsData$LogWage2)

# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL           13.1762012       0                   1     1.37391920
# KS            1.4650379       0                   1     0.68258658
# CvM           0.7557939       0                   1     0.06737924
# AD           23.4758645       0                   1     3.16250502
# spearman      0.5585395       0                   1     0.13510847
# kendall       0.4210976       0                   1     0.09105750
