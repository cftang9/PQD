#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
###########################################################################
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
EducationData = read.csv("https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.csv")
# sample size n=51; 
n = length(EducationData$GraRate); 

###########################################################################
# Here we provide the scatterplot of the data and pseudo-observations 
# to roughly visualize the dependence structure between 
# graduation rate and amount spent per person. 
par(mar=c(4.5,5,3,0.5))
par(mfrow=c(1,2))
plot(EducationData$GraRate, EducationData$SpentStud,xlab="Graduation rate",ylab="Amount spent per student",main="Scatterplot of the data") 
plot(rank(EducationData$GraRate)/(n+1),rank(EducationData$SpentStud)/(n+1),xlab="Graduation rate",ylab="Amount spent per student",main="Scatterplot of pseudo-observations"); 

# Perform all considered tests: 
# It takes less than 2 mins
set.seed(100)
IndvsPQD(EducationData$GraRate,EducationData$SpentStud)

# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL           1.14408178  0.0958                   0     1.43299162
# KS           0.47067846  0.3301                   0     0.67176676
# CvM          0.08260912  0.0678                   0     0.09068779
# AD           2.77333756  0.2244                   0     4.71607090
# spearman     0.21485303  0.0682                   0     0.23475566
# kendall      0.15015742  0.0615                   0     0.15921569
