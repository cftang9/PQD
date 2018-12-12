#########################################################################
#   R functions for testing independence versus positive quadrant       #
# dependence corresponding to the manuscript titled,                    #
#             "Testing for positive quadrant dependence."               #
#                         Date: 06/17/2018                              #
#########################################################################
library(Rcpp)
library(copula)
source("http://people.stat.sc.edu/wang528/PQD/EL_PQD_Library.R")
##########################################################################
#                                                                        #
#    For illustration, we generate a random sample from a Clayton copula #
#    with a user-specifiic Kandell's tau test for independence versus    #
#    positive quadrant dependence.                                       #
#    (For large sample sizes, it may take some time to generate the      #
#     critical values)                                                   #
#                                                                        #
##########################################################################

n = 10; X = array(,c(n,2)); 
tau = 0.2; 
set.seed(100);
X = RV_CopTau(n, tau, Copula="Clayton"); 
IndvsPQD(X=X[,1],Y=X[,2],graph=TRUE)
# [1] "1: reject independence; 0: do not rejct independence"
#           test_statistic p-value reject_independence critical_value
# EL           0.39887816  0.5200                   0      1.4329523
# KS           0.31884122  0.8956                   0      0.6664304
# CvM          0.03267605  0.8528                   0      0.1961564
# AD           1.98834920  0.7074                   0      7.8084519
# spearman    -0.17575758  0.6902                   0      0.5515152
# kendall     -0.20000000  0.7611                   0      0.4222222


##########################################################################
#                                                                        #
#    For illustration, we generate a random sample from a Frank copula   #
#    with a user-specifiic Kandell's tau test for independence versus    #
#    positive quadrant dependence.                                       #
#    (For large sample sizes, it may take some time to generate the      #
#     critical values)                                                   #
#                                                                        #
##########################################################################

n = 10; X = array(,c(n,2)); 
tau = 0.4; 
set.seed(100);
X = RV_CopTau(n, tau, Copula="Frank"); 
IndvsPQD(X=X[,1],Y=X[,2],graph=TRUE)
# [1] "1: reject independence; 0: do not rejct independence"
#           test_statistic p-value reject_independence critical_value
# EL            2.9119301  0.0001                   1      1.4329523
# KS            0.8572125  0.0026                   1      0.6664304
# CvM           0.4161669  0.0000                   1      0.1961564
# AD           13.8342732  0.0002                   1      7.8084519
# spearman      0.9151515  0.0001                   1      0.5515152
# kendall       0.7777778  0.0000                   1      0.4222222

##########################################################################
#                                                                        #
#    For illustration, we generate a random sample from a Gumbel copula  #
#    with a user-specifiic Kandell's tau test for independence versus    #
#    positive quadrant dependence.                                       #
#    (For large sample sizes, it may take some time to generate the      #
#     critical values)                                                   #
#                                                                        #
##########################################################################

n = 10; X = array(,c(n,2)); 
tau = 0.4; 
set.seed(100);
X = RV_CopTau(n, tau, Copula="Gumbel"); 
IndvsPQD(X=X[,1],Y=X[,2],graph=TRUE)
# [1] "1: reject independence; 0: do not rejct independence"
#          test_statistic p-value reject_independence critical_value
# EL            1.1198916  0.1065                   0      1.4629308
# KS            0.5331443  0.3095                   0      0.6664304
# CvM           0.1284072  0.1957                   0      0.1985585
# AD            4.8414675  0.2288                   0      7.7588699
# spearman      0.4545455  0.0882                   0      0.5515152
# kendall       0.2888889  0.1058                   0      0.4222222

##########################################################################
#                                                                        #
#    For illustration, we generate a random sample from a bivariate      #
#    Gaussian copula with a user-specifiic correlation rho to            #
#    test for independence versus positive quadrant dependence.          #
#    (For large sample sizes, it may take some time to generate the      #
#     critical values)                                                   #
#                                                                        #
##########################################################################

n = 10; X = array(,c(n,2));
rho = 0.0; 
set.seed(100);
X = RV_CopGaussian(n, rho)
IndvsPQD(X=X[,1],Y=X[,2],graph=TRUE)
# [1] "1: reject independence; 0: do not rejct independence"
#           test_statistic p-value reject_independence critical_value
# EL           0.15735832  0.8154                   0      1.4329523
# KS           0.32668158  0.8662                   0      0.6664304
# CvM          0.03423537  0.8401                   0      0.1961564
# AD           1.27527351  0.8626                   0      7.8084519
# spearman    -0.28484848  0.7952                   0      0.5515152
# kendall     -0.20000000  0.7610                   0      0.4222222