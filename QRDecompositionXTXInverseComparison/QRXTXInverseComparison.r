# Compare efficiency of three QR decomposition related methods of computing inverse(X'X)
# Author:  Tom Balmat
# Date:    Nov 2017

# Let X be a (n x p) design matrix
# Let Qx and Rx be the QR decomposition of X
# Let Qxtx and Rxtx be the QR decomposition of X'X

# Note the following properties of a QR decomposition:
# 1.) Q'Q = I
# 2.) R is upper triangular (and, therefore, assumed invertible)
# 3.) QR = X => R'Q'QR = R'R = X'X
# 4.) inverse(X'X)*R'R = I => inverse(R) = inverse(X'X)*R' =>
#     inverse(X'X) = inverse(R)*inverse(R')
# 5.) inverse(R') = [inverse(R)]' => inverse(X'X) = inverse(R)*[inverse(R)]'
# 6.) X'X = R'R => inverse(X'X)=inverse(R'R)

# Compare time to execute the following
# 1.) compute QR(X), compute inverse(X'X) = inverse(R'R)
# 2.) compute QR(X), compute inverse(X'X) = inverse(R)*[inverse(R)]'
# 3.) compute X'X, compute QR(X'X), compute inverse(X'X) from solve(QR)

options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)

setwd("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeModelSolutionComparison")

#####################################################################################################
#### load Matrix libraries for efficient X'X construction
#####################################################################################################

library(Matrix)
library(MatrixModels)

#####################################################################################################
#### function to generate fixed effects test data
#####################################################################################################

feTestData <- function(n, nlv3, nlv4, nlv5) {
 
  # model:  Y = b0 + b1X1 + b2X2 + b3iX3i + b4jX4j + b5kX5k + e~N(0, evar)

  # function parameters:
  # n ..... number of observations to generate
  # nlv3 .. number of levels to generate for fixed effect 3
  # nlv4 .. number of levels to generate for fixed effect 4
  # nlv5 .. number of levels to generate for fixed effect 5

  # static parameters (modify as needed)
  # b0 .... intercept
  # b1 .... coefficient of X1
  # X1 .... vector of equally weighted values to sample for X1
  # b2 .... coefficient of X2
  # X2 .... vector of equally weighted values to sample for X2
  # b3 .... vector of X3 fixed effect coefficients (length(b3) equally weighted levels will be generated)
  # b4 .... vector of X4 fixed effect coefficients (length(b4) equally weighted levels will be generated)
  # b5 .... number of X5 fixed effect coefficients (length(b5) equally weighted levels will be generated)
  # sd .... standard deviation of normally distributed errors
  b0 <- 10.333
  b1 <- 2.3
  X1 <- sample(seq(2, 4, 0.01), n, replace=T)
  b2 <- 5.5
  X2 <- sample(seq(3, 5, 0.01), n, replace=T)
  b3 <- sample(seq(0.5, 2.5, 2/nlv3)[1:nlv3])
  X3 <- sample(1:length(b3), n, replace=T)
  b4 <- sample(seq(2.5, 5, 2.5/nlv4)[1:nlv4])
  X4 <- sample(1:length(b4), n, replace=T)
  b5 <- sample(seq(5, 7.5, 2.5/nlv5)[1:nlv5])
  X5 <- sample(1:length(b5), n, replace=T)
  sd <- 10

  # set reference level coefficients to 0 so that b0 and ref levels unconfounded
  # this enables comparison of estimate b0 to known b0
  b3[1] <- 0
  b4[1] <- 0
  b5[1] <- 0

  # generate simulated indepent and calculated dependent values
  data.frame("Y" = rep(b0, n) + b1*X1 + b2*X2 + b3[X3] + b4[X4] + b5[X5] + rnorm(n, 0, sd),
             "X1"=X1, "X2"=X2,
             "X3"=sprintf("%05.0f", X3), "X4"=sprintf("%05.0f", X4), "X5"=sprintf("%05.0f", X5))
       
}

#####################################################################################################
#### configure test parameters
#####################################################################################################

dpar <- rbind(data.frame("n"=10000,    "lv3"=10,   "lv4"=15,   "lv5"=20),
              data.frame("n"=25000,    "lv3"=10,   "lv4"=20,   "lv5"=40),
              data.frame("n"=50000,    "lv3"=15,   "lv4"=30,   "lv5"=50),
              data.frame("n"=100000,   "lv3"=15,   "lv4"=40,   "lv5"=60),
              data.frame("n"=150000,   "lv3"=25,   "lv4"=40,   "lv5"=60),
              data.frame("n"=250000,   "lv3"=25,   "lv4"=50,   "lv5"=75),
              data.frame("n"=500000,   "lv3"=50,   "lv4"=75,   "lv5"=100),
              data.frame("n"=1000000,  "lv3"=75,   "lv4"=150,  "lv5"=250))
#              data.frame("n"=1500000,  "lv3"=100,  "lv4"=250,  "lv5"=500),
#              data.frame("n"=2000000,  "lv3"=100,  "lv4"=500,  "lv5"=1000),
#              data.frame("n"=5000000,  "lv3"=200,  "lv4"=1000, "lv5"=5000),
#              data.frame("n"=10000000, "lv3"=1000, "lv4"=5000, "lv5"=10000),
#              data.frame("n"=25000000, "lv3"=1000, "lv4"=5000, "lv5"=25000))
              
#####################################################################################################
#### execute five iterations of each method for each parameter set
#####################################################################################################

for(i in 1:5)
  for(j in 1:nrow(dpar)) {
    n <- dpar[j,"n"]
    lv3 <- dpar[j,"lv3"]
    lv4 <- dpar[j,"lv4"]
    lv5 <- dpar[j,"lv5"]

    # generate sample fixed effects data set
    feDat <- feTestData(n, lv3, lv4, lv5)

    # generate X and X'X by expanding fixed effects into indicator columns, one for each
    # non-reference level of each fixed effect
    X0 <- sparse.model.matrix(~.,
            data=data.frame(
              feDat[,c("X1", "X2")],
              "X3"=relevel(factor(feDat[,"X3"]), ref=min(feDat[,"X3"])),
              "X4"=relevel(factor(feDat[,"X4"]), ref=min(feDat[,"X4"])),
              "X5"=relevel(factor(feDat[,"X5"]), ref=min(feDat[,"X5"])) ))
    X <- as.matrix(X0)
    rm(feDat)
    gc()

    # compute QR decomposition of X
    t <- proc.time()[["elapsed"]]
    qrX <- qr(X)
    tqr <- proc.time()[["elapsed"]]-t
    gc()

    # compute inverse(X'X) using inverse(R)*[inverse(R)]'
    t <- proc.time()[["elapsed"]]
    iR <- solve(qr.R(qrX))
    iXTX <- iR%*%t(iR)
    qrd <- data.frame("method"="iRt[iR]", "p"=ncol(iXTX), "t"=tqr+proc.time()[["elapsed"]]-t)
    write.table(qrd, "QRXTXInverseComparison.csv", sep=", ", row.names=F, col.names=F, quote=F, append=T)
    rm(list=c("iR", "iXTX", "qrd"))
    gc()

    # compute inverse(X'X) using inverse(R'R)
    t <- proc.time()[["elapsed"]]
    R <- qr.R(qrX)
    iXTX <- solve(t(R)%*%R)
    qrd <- data.frame("method"="iRTR", "p"=ncol(iXTX), "t"=tqr+proc.time()[["elapsed"]]-t)
    write.table(qrd, "QRXTXInverseComparison.csv", sep=", ", row.names=F, col.names=F, quote=F, append=T)
    rm(list=c("R", "iXTX", "qrd"))
    gc()

    # compute X'X, QR decomposition of X'X, then inverse(X'X) using qr.solve(QR decomp)
    # note the masking of %*% to efficient sparse matrix operations
    XTX <- as.matrix(t(X0)%*%X0)
    t <- proc.time()[["elapsed"]]
    qrXTX <- qr(XTX)
    iXTX <- qr.solve(qrXTX)
    qrd <- data.frame("method"="qrXTX", "p"=ncol(iXTX), "t"=proc.time()[["elapsed"]]-t)
    write.table(qrd, "QRXTXInverseComparison.csv", sep=", ", row.names=F, col.names=F, quote=F, append=T)
    rm(list=c("XTX", "qrXTX", "iXTX", "qrd"))
    gc()

  }

#which(abs(iXTX0-iXTX)>1e15)
