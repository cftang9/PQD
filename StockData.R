#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
source("https://raw.githubusercontent.com/cftang9/PQD/master/EL_PQD_Library.R")
########################################################################
#   This data is download from https://finance.yahoo.com/ 
StockData = read.csv("https://raw.githubusercontent.com/cftang9/PQD/master/StockData.csv",head=TRUE)[,2:4]
n = length(StockData$AAPL); #sample size n=124
##########################################################################
# Perform Ljung–Box test to examine the independence of the returns.
# the resulting p-values are all larger than alpha=0.05.

Box.test( StockData$AAPL , type='Ljung')
# Box-Ljung test
# data:  StockData$AAPL
# X-squared = 0.31136, df = 1, p-value = 0.5768

Box.test( StockData$GOOGL , type='Ljung')
# Box-Ljung test
# data:  StockData$GOOGL
# X-squared = 0.22509, df = 1, p-value = 0.6352

Box.test( StockData$WMT , type='Ljung')
# Box-Ljung test
# data:  StockData$WMT
# X-squared = 3.0171, df = 1, p-value = 0.08239


###########################################################################
# Here we provide the Scatterplot of the data and pseudo-observations 
# to roughly visualize the dependence structure between 
# each pair of the three stocks. 
par(mfcol=c(2,3))
par(mar=c(4.5,5,3,0.5))
plot(StockData$AAPL, StockData$GOOGL, xlab="GOOGL", ylab="AAPL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of the data")
plot(rank(StockData$AAPL)/(n+1), rank(StockData$GOOGL)/(n+1), xlab="GOOGL", ylab="AAPL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of pseudo-observations")

plot(StockData$AAPL, StockData$WMT, xlab="WMT", ylab="AAPL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of the data")
plot(rank(StockData$AAPL)/(n+1), rank(StockData$WMT)/(n+1), xlab="WMT", ylab="AAPL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of pseudo-observations")

plot(StockData$GOOGL, StockData$WMT, xlab="WMT", ylab="GOOGL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of the data")
plot(rank(StockData$GOOGL)/(n+1), rank(StockData$WMT)/(n+1), xlab="WMT", ylab="GOOGL", cex.lab=1.25, cex.axis=1.25,main="Scatterplot of pseudo-observations")

# Perform all considered tests: 

# AAPL VS GOOGL
# It takes <= 10 mins
set.seed(100)
IndvsPQD(StockData$AAPL,StockData$GOOGL)

# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL            8.2532935   0e+00                   1     1.40294861
# KS            1.0180776   7e-04                   1     0.67873562
# CvM           0.3808337   0e+00                   1     0.06943679
# AD           17.4198196   0e+00                   1     3.36676649
# spearman      0.4762549   0e+00                   1     0.14824579
# kendall       0.3406766   0e+00                   1     0.09967217

  
# AAPL VS WMT
# It takes <= 10 mins
set.seed(100)
IndvsPQD(StockData$AAPL,StockData$WMT)

# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL           1.40990425  0.0495                   1     1.40294861
# KS           0.66164869  0.0587                   0     0.67873562
# CvM          0.08895359  0.0213                   1     0.06943679
# AD           3.80117839  0.0317                   1     3.36676649
# spearman     0.15488277  0.0418                   1     0.14824579
# kendall      0.11224757  0.0312                   1     0.09967217

# GOOGL VS WMT
# It takes <= 10 mins
set.seed(100)
IndvsPQD(StockData$GOOGL,StockData$WMT)

# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL            4.2656097  0.0001                   1     1.40294861
# KS            1.0279458  0.0004                   1     0.67873562
# CvM           0.2465081  0.0000                   1     0.06943679
# AD            6.8627437  0.0016                   1     3.36676649
# spearman      0.3359371  0.0001                   1     0.14824579
# kendall       0.2355101  0.0001                   1     0.09967217
