# Duke University SSRI Human Capital Project
# Large Fixed Effects Model Solution
# Comparison of fixed effects solution algorithms using OPM data set

options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)

library(RODBC)             # database i/f

# specify current working dir (fixed effect regression source script location)
setwd("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeModelSolutionComparison\\")

#####################################################################################################
#### compare lfe, feXTX, and sparse matrix functions from Matrix and Matrix Models packages
#####################################################################################################

# query source data
db <- odbcDriverConnect(connection="driver={SQL Server}; server=ssri-swrk-pap36; database=HC_Dev; trusted_connection=true", readOnlyOptimize=T)
t <- proc.time()
opm <- sqlQuery(db, "FixedEffectsModelData @style='NonUniqueIDJdF', @WorkSchedule='FullTime'", stringsAsFactors=F)
proc.time()-t
odbcClose(db)

####
#### method:  lfe()
####

feSolution_lfe <- function(opmDat) {
  # fit linear fixed effects model using lfe()
  # parameters:
  # opmDat .... source OPM data (lnBasicPay, age, age**2, education, sex, race, year, bureau, occupation)
  library(lfe)
  # convert fixed effects to factors
  # note that fixed effects involved in interaction terms must be factors
  opmDat <- data.frame(opmDat[,c("lnBasicPay", "Age", "AgeSq", "EducationYears")],
                       "Sex"=relevel(factor(opmDat[,"Sex"]), ref="M"),
                       "Race"=relevel(factor(opmDat[,"Race"]), ref="E"),
                       "FY"=relevel(factor(opmDat[,"FY"]), ref="1988"),
                       "BureauID"=relevel(factor(opmDat[,"BureauID"]), ref="114009000"),
                       "Occupation"=relevel(factor(opmDat[,"Occupation"]), ref="303"))
  # fit model
  t <- proc.time()
  m <- felm(lnBasicPay ~ Age+AgeSq+EducationYears | Sex+Race+Sex:Race+FY+BureauID+Occupation, data=opmDat, exactDOF='rM')
  # following seems to ignore method='cg' and imposes default Kaczmarz method
  # reference levels are ignored and first level of FE 1 is coerced into a reference level for all FEs
  # also, an intercept term is not estimated
  b <- getfe(m, se=T, ef=efactory(m))
  # return parameter estimates and compute time
  # note the absence of an intercept term
  list("est"=c("b0"=0,
               m$coefficients[c("Age", "AgeSq", "EducationYears"), "lnBasicPay"],
               b[c("Sex.F", "Race.A", "Race.B", "Race.C", "Race.D"), "effect"],
               b[c("Sex:Race.F.A", "Sex:Race.F.B", "Sex:Race.F.C", "Sex:Race.F.D"), "effect"],
               b[which(b[,"fe"]=="FY" & b[,"idx"]!="1988"),"effect"],
               b[which(b[,"fe"]=="BureauID" & b[,"idx"]!="114009000"),"effect"],
               b[which(b[,"fe"]=="Occupation" & b[,"idx"]!="303"),"effect"]),
       "time"=(proc.time()-t)[1:3])
}

####
#### method:  feXTX
####

feSolution_feXTX <- function(opmDat, ncore=12) {
  # fit linear fixed effects model using sparse XTX indicator matrix construction
  # parameters:
  # opmDat .... source OPM data (lnBasicPay, age, age**2, education, sex, race, year, bureau, occupation)
  # ncores .... number of parallel cores to use in composing sparse matrix
  # load Cholesky and feXTX functions
  sys.source("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\FixedEffectsMatrixSolution.r", env=environment())
  sys.source("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\choleskyDecompLoad.r", env=environment())
  sys.source("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\cholInvDiagLoad.r", env=environment())
  t <- proc.time()
  # compose interaction indicators prior to feXTX() call
  fXraceA <- rep(0, nrow(opmDat))
  fXraceA[which(opmDat[,"Sex"]=="F" & opmDat[,"Race"]=="A")] <- 1
  fXraceB <- rep(0, nrow(opmDat))
  fXraceB[which(opmDat[,"Sex"]=="F" & opmDat[,"Race"]=="B")] <- 1
  fXraceC <- rep(0, nrow(opm))
  fXraceC[which(opmDat[,"Sex"]=="F" & opmDat[,"Race"]=="C")] <- 1
  fXraceD <- rep(0, nrow(opm))
  fXraceD[which(opmDat[,"Sex"]=="F" & opmDat[,"Race"]=="D")] <- 1
  # fit model
  # note the specification of Cholesky decomposition parallelized methods beta and variance estimation
  m <- feXTX(data=data.frame(opmDat, "fXraceA"=fXraceA, "fXraceB"=fXraceB, "fXraceC"=fXraceC, "fXraceD"=fXraceD),
             Y="lnBasicPay",
             contX=c("Age", "AgeSq", "EducationYears", "fXraceA", "fXraceB", "fXraceC", "fXraceD"),
             fixedX=c("Sex", "Race", "FY", "BureauID", "Occupation"),
             refLevel=c("M", "E", "1988", "114009000", "303"),
             interactionX=NULL,
             estBetaVar="stdOLS",
             robustVarID=NULL,
             nCoreXTX=ncore,
             nCoreVar=0,
             solMethod="chol-parallel")
  # return list of estimates (in standard order) and execution time
  list("est"=c(m$beta[c("b0", "x-Age", "x-AgeSq", "x-EducationYears",
                        "Sex-F", "Race-A", "Race-B", "Race-C", "Race-D",
                        "x-fXraceA", "x-fXraceB", "x-fXraceC", "x-fXraceD")],
               m$beta[which(substring(names(m$beta), 1, 2)=="FY")],
               m$beta[which(substring(names(m$beta), 1, 8)=="BureauID")],
               m$beta[which(substring(names(m$beta), 1, 10)=="Occupation")]),
       "time"=(proc.time()-t)[1:3])
}

####
#### method:  Matrix package, sparse matrix, lm.fit and solve from package
####

feSolution_MatrixSparseFitSolve <- function(opmDat) {
  # fit linear fixed effects model using sparse matrix from Matrix package, solution using lm.fit from package
  # parameter variances estimated using solve function from package
  # parameters:
  # opmDat .... source OPM data (lnBasicPay, age, age**2, education, sex, race, year, bureau, occupation)
  library(Matrix)
  library(MatrixModels)
  t <- proc.time()
  # construct sparse matrix
  X <- sparse.model.matrix(~Age+AgeSq+EducationYears+FY+Sex+Race+Sex:Race+BureauID+Occupation,
                           data=data.frame(opmDat[,c("Age", "AgeSq", "EducationYears")],
                                           "Sex"=relevel(factor(opmDat[,"Sex"]), ref="M"),
                                           "Race"=relevel(factor(opmDat[,"Race"]), ref="E"),
                                           "FY"=relevel(factor(opmDat[,"FY"]), ref="1988"),
                                           "BureauID"=relevel(factor(opmDat[,"BureauID"]), ref="114009000"),
                                           "Occupation"=relevel(factor(opmDat[,"Occupation"]), ref="303")))
  # solve parameter estimates
  beta <- MatrixModels:::lm.fit.sparse(X, opmDat[,"lnBasicPay"], method="cholesky")$coef
  names(beta) <- colnames(X)
  # compute estimate of Y variance
  yVar <- sum((opmDat[,"lnBasicPay"]-X%*%beta)**2)/(nrow(X)-ncol(X))
  # compute estimate of beta variance
  vBeta <- yVar*diag(Matrix::solve(Matrix::crossprod(X)))
  # return list of estimates and execution time
  list("est"=c(beta[1],
               beta[c("Age", "AgeSq", "EducationYears")],
               beta[c("SexF", "RaceA", "RaceB", "RaceC", "RaceD")],
               beta[c("SexF:RaceA", "SexF:RaceB", "SexF:RaceC", "SexF:RaceD")],
               beta[which(substring(names(beta), 1, 2)=="FY")],
               beta[which(substring(names(beta), 1, 8)=="BureauID")],
               beta[which(substring(names(beta), 1, 10)=="Occupation")]),
       "time"=(proc.time()-t)[1:3])
}

####
#### method:  Matrix package, sparse matrix, custom parallel Cholesky functions
####

feSolution_MatrixSparseChol <- function(opmDat) {
  # fit linear fixed effects model using sparse matrix from Matrix package, solution using custom parallel Cholesky
  # algorithm, parameter variances estimated using custom parallel inverse from Cholesky decomposition algorithm
  # parameters:
  # opmDat .... source OPM data (lnBasicPay, age, age**2, education, sex, race, year, bureau, occupation)
  library(Matrix)
  library(MatrixModels)
  # load Cholesky functions
  sys.source("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\choleskyDecompLoad.r", env=environment())
  sys.source("\\\\SSRI-NAS-FE01.oit.duke.edu\\SSRI\\OPM\\Users\\Current\\tjb48\\Analysis\\FixedEffectsModel\\LargeFixedEffectsModel\\CholeskyDecomposition\\cholInvDiagLoad.r", env=environment())
  t <- proc.time()
  # construct sparse matrix
  X <- sparse.model.matrix(~Age+AgeSq+EducationYears+FY+Sex+Race+Sex:Race+BureauID+Occupation,
                           data=data.frame(opmDat[,c("Age", "AgeSq", "EducationYears")],
                                           "Sex"=relevel(factor(opmDat[,"Sex"]), ref="M"),
                                           "Race"=relevel(factor(opmDat[,"Race"]), ref="E"),
                                           "FY"=relevel(factor(opmDat[,"FY"]), ref="1988"),
                                           "BureauID"=relevel(factor(opmDat[,"BureauID"]), ref="114009000"),
                                           "Occupation"=relevel(factor(opmDat[,"Occupation"]), ref="303")))
  # compute Cholesky decomposition of X'X
  chXTX <- choleskyDecomp(as.matrix(t(X)%*%X))
  # compute X'Y
  XTY <- as.vector(t(X)%*%opmDat[,"lnBasicPay"])
  # solve parameter estimates
  beta <- backsolve(chXTX, forwardsolve(chXTX, XTY, upper.tri=T, transpose=T), upper.tri=T)
  names(beta) <- colnames(X)
  # compute estimate of Y varaiance
  yVar <- sum((opmDat[,"lnBasicPay"]-X%*%beta)**2)/(nrow(X)-ncol(X))
  # compute estimate of beta variance
  vBeta <- yVar*diag(cholInvDiag(chXTX))
  # return list of estimates and execution time
  # 
  list("est"=c(beta[1],
               beta[c("Age", "AgeSq", "EducationYears")],
               beta[c("SexF", "RaceA", "RaceB", "RaceC", "RaceD")],
               beta[c("SexF:RaceA", "SexF:RaceB", "SexF:RaceC", "SexF:RaceD")],
               beta[which(substring(names(beta), 1, 2)=="FY")],
               beta[which(substring(names(beta), 1, 8)=="BureauID")],
               beta[which(substring(names(beta), 1, 10)=="Occupation")]),
       "time"=(proc.time()-t)[1:3])
}

####
#### execute methods and retain results
####

sollfe <- feSolution_lfe(opm)
gc()
solfeXTX <- feSolution_feXTX(opm, ncore=20)
gc()
solMatrix1 <- feSolution_MatrixSparseFitSolve(opm)
gc()
solMatrix2 <- feSolution_MatrixSparseChol(opm)
gc()

####
#### inspect compute time in minutes
####

sollfe$time/60
solfeXTX$time/60
solMatrix1$time/60
solMatrix2$time/60

####
#### compose table of coefficient estimate differences (to those of feXTX)
####

bdiff <- t(apply(as.matrix(cbind(sollfe$est, solMatrix1$est, solMatrix2$est)), 2,
                 function(b) {
                   d <- abs(b-solfeXTX$est)
                   c(length(which(d<1e-12)),
                     length(which(d>=1e-12 & d<1e-10)),
                     length(which(d>=1e-10 & d<1e-7)),
                     length(which(d>=1e-7 & d<1e-5)),
                     length(which(d>=1e-5)))
                 }))
rownames(bdiff) <- c("lfe", "Matrix1", "Matrix2")
colnames(bdiff) <- c("0<1e-12", "1e-12<1e-10", "1e-10<1e-7", "1e-7<1e-5", ">1e-5")


