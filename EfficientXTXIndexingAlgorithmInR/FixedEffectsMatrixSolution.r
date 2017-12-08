# Duke University Human Capital Project
# Duke University Synthetic Data Project
# Efficient X'X construction Method for Solution of Large Observation, High Dimension Fixed Effects Models

# Author:  Tom Balmat
# Version:  5/14/2017

# Copyright 2017 Duke University, Durham, North Carolina

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
# following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial
# portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

feXTX <- function(data, Y, contX, fixedX, refLevel, interactionX=NULL, estBetaVar="none", robustVarID=NULL,
                  nCoreXTX, nCoreVar, solMethod="chol") {

  ##############################################################################################################
  #### construct X'X and X'Y using variable cross products using sparse vector indices
  #### solve X'X*Beta = X'Y for Beta using qr decomposition
  #### estimate Beta standard errors using classic OLS inv(X'X), robust, or clustered based on specified ID 
  ####
  #### fixed effect factors are expanded into a list of vectors, one for each non-reference level of the factor,
  #### each consisting of pointers to observations corresponding to the associated factor level
  ####
  #### factor level pointers are used in simple sum operations as an efficient alternative to exhaustive
  #### row-transpose, column products and is very efficient in high degree low density factor situations
  ####
  #### library dependencies:  parallel (has been successfully tested with SNOW, also)
  ####
  #### this implementation is designed for MS Windows, which does not offer shared memory support, and involves
  #### (expensive execution and memory) export of objects to individual core memory for parallel operations
  ####
  #### for Unix implementation:
  #### create a "forked" type socket cluster, makeCluster(type="FORK"), to enable shared access to memory and
  #### omit clusterExport operations that create unnecessary copies of large vectors
  ####
  #### parameters:
  #### data:           data frame containing continuous, fixed effect, and dependent variables
  ####                 note that if fixed effects are supplied as factors then levels become factor values
  ####                 these are easily translated to original values with X[,levels(X)[i]]
  #### Y:              name of column in data containing dependent numeric values
  #### contX:          vector of column names of numeric data containing continuous independent values
  ####                 specify contX=NULL if  there are no continuous variables in the model
  #### fixedX:         vector of column names of data containing fixed effect independent values
  ####                 must be non-NULL (this is a fixed effects solution!)
  #### refLevel:       vector of fixed effect reference levels, in order of fixedX columns
  #### interactionX:   two column data frame or matrix of interaction pairs (each row constitutes a pair)
  ####                 note that if an interaction involves a continuous variable and a fixed effect then
  ####                 the continuous variable is moved to the left side if not already declared on the left
  ####                 this reduces the number of possible continuous-fixed interaction pairings from four
  ####                 (cc, cf, fc, ff) to three (cc, cf, ff)
  ####                 specify interactionX=NULL if  there are no interactions in the model
  #### estBetaVar:     beta variance estimation method
  ####                   "stdOLS"    for standard OLS estimates [diagonal of inv(X'X)]
  ####                   "robust"    for robust variance (independent, non-identically-distributed errors)
  ####                   "cluster"   for variance based on clustering around IDs specified in robustVarID
  ####                   "i(x'x)x'uu'xi(x'x)" to confirm analytical result of zero (see notes under variances)
  ####                   "none"      for no variance estimates
  ####                 apologies, but robust, clustered, and i(x'x)x'uu'xi(x'x) variances are computed only
  ####                 for models with no interactions specified (this feature will be added in time)
  #### robustVarID:    name of column in data containing vector of IDs used in clustered variance estimatation
  #### nCoreXTX:       the number of cores to use in parallel XTX composition operations
  #### nCoreVar:       the number of cores to use in parallel robust/clustered variance estimation operations
  #### solMethod:      method of solving X'Xb=X'Y
  ####                   "qr" for QR decomposition
  ####                   "chol" for Cholesky decomposition (using R chol() function)
  ####                   "chol-parallel" for parallel Cholesky decomposition (using custom choleskyDecomp
  ####                     function) - this also coerces stdOLS parameter variance computation using custom
  ####                     parallel cholInvDiag function 
  ####
  #### the return object is a list with the following elements:
  #### status:          completion status ("" or error message)
  #### Y:               echo of input parameter
  #### contX:           echo of input parameter
  #### fixedX:          echo of input parameter
  #### refLevel:        echo of input parameter
  #### interactionX:    echo of input parameter
  #### nX:              number of rows in data
  #### beta:            OLS parameter estimates for Y ~ b0 + contX + indicator(fixedX|not ref)
  #### vBeta:           parameter estimate variances
  #### estBetaVar:      echo of input parameter
  #### XTX:             X'X matrix (the primary effort of this function!)
  #### dimXTX:          column (row) dimension of X'X
  #### rankXTX:         rank of X'X as reported by qr() or chol()
  #### XTY:             X'Y
  #### yEst:            estimated Y values
  #### yVar:            variance of observed Y's
  #### time:            matrix of execution time, row 1 for compute, row 2 for memory export
  ##########################################################################################################

  ##########################################################################################################
  # to be implemented:
  # verify that specified fixed effect reference levels exists in supplied data
  # otherwise, a linear dependent condition exists between fixed effect(s) and the constant column
  # verify that interaction variables appear in data
  # verify that 0 < nCoreXTX and nCoreVar < available cores
  # convert order() calls in fixed effect indicator expansion to sort.list(factor(x), method="quick", na.last=NA)
  ##########################################################################################################

  library(parallel)

  ###########################################################################################
  #### configure program objects and parameters
  ###########################################################################################

  tcompute <- rep(0, 5)
  tmem <- rep(0, 5)
  t0 <- proc.time()

  # verify that specified dependent, continuous, and fixed effect names are valid (length 1)
  # verify that specified dependent, continuous, and fixed effect vectors exist in supplied data
  # dependent var
  if(class(Y)=="character") {
    if(length(Y)==1) {
      if(length(which(colnames(data)==Y))!=0) {
        status <- ""
      } else {
        status <- paste("ERROR (feXTX) - Specified Y vector missing in source data (", Y, ")", sep="")
      }
    } else {
      status <- "ERROR (feXTX) - Invalid Y specification (single character string required)"
    }
  } else {
    status <- "ERROR (feXTX) - Invalid Y specification (single element character vector required)"
  }
  # continuous vars
  # note that continuous vars are not a requirement (models may be specified with none)
  if(status=="" & !is.null(contX)) {
    if(class(contX)=="character") {
      if(length(contX)>0) {
        cid <- which(!contX %in% colnames(data))
        if(length(cid)==0) {
          status <- ""
        } else {
          status <- paste("ERROR (feXTX) - Specified contX vector(s) missing in source data (",
                          paste(contX[cid], collapse=" "), ")", sep="")
        }
      }
    } else {
      status <- "ERROR (feXTX) - Invalid contX specification (vector of character names or NULL required)"
    }
  }
  # fixed effects
  if(status=="") {
    if(class(fixedX)=="character") {
      if(length(fixedX)>0) {
        cid <- which(!fixedX %in% colnames(data))
        if(length(cid)==0) {
          status <- ""
        } else {
          status <- paste("ERROR (feXTX) - Specified fixedX vector(s) missing in source data (",
                          paste(fixedX[cid], collapse=" "), ")", sep="")
        }
      } else {
        status <- "ERROR (feXTX) - Invalid fixedX specification (at least one fixed effect must be specified)"
      }
    } else {
      status <- "ERROR (feXTX) - Invalid fixedX specification (vector of character names required)"
    }
  }
  
  if(status=="") {
  
    # count levels of fixed effects for ordering
    # moving high dimension effects to the end of list improves performance in paralell apply
    # operations since this distributes more levels to available cores in a give cross multiply
    fixedOrder <- order(apply(as.matrix(1:length(fixedX)), 1, function(j) length(unique(data[,fixedX[j]]))))
    fixedX <- fixedX[fixedOrder]
    refLevel <- refLevel[fixedOrder]

    # save input parameters (they are reported in the function output)
    fparams <- list("Y"=Y, "contX"=contX, "fixedX"=fixedX)

    # save row dimension and number of continuous cols of design matrix
    # note that NULL contX has length of 0
    nX <- nrow(data)
    ncontX <- length(contX)
    
    # convert independent and dependent vars to numeric to avoid potential integer overflow
    if(length(contX)>0)
      contX <- apply(as.matrix(data[,contX]), 2, as.numeric)
    fixedX <- data.frame(data[,fixedX])
    colnames(fixedX) <- fparams[["fixedX"]]
    Y <- as.numeric(data[,Y])

    # continue only if fixed effects declared
    if(ncol(fixedX)>0) {

      ##############################################################################################
      #### expand fixed effects into indicator columns for efficient (sparse, non-zero) addressing
      ##############################################################################################

      # construct sorted unique fixed effect value lists, excluding specified reference levels
      uFixedLevel <- lapply(1:ncol(fixedX), FUN=function(i) {
                                                  x <- sort(unique(fixedX[,i]))
                                                  x[which(x!=refLevel[i])]
                                                })

      # count levels by fixed effect
      nFixedLevel <- mapply(1:length(uFixedLevel), FUN=function(i) length(uFixedLevel[[i]]))

      # define function to execute parallel operations
      # this enables assigning the global environment as its env, which avoids parLapply's
      # insistence on exporting the entire local environment (it never exports .GlobalEnv,
      # why it exports anything is a puzzle since environment objects cannot be referenced
      # in parallel functions anway) 
      # this is important if executing within a function with a local environment
      # since it prevents export of potentially large objects
      #plyXTX <- function(cl, ulev){
        # unique FE levels are distributed to cores, each sub-process returns a
        # vector of row positions that contain its corresponding level
        # parameters:
        # cl ............ cluster of cores on which to distribute individual lapply instances
        # uFixedLevel ... list of unique levels to identify observations by (using which)
        #                 levels are sequentially passed to the function as lev
        # function ...... apply which to fixed effects column searching for lev and return vector
        #                 of observation indices that contain lev
        # the return object is a list of which() results (vectors), one for each level in sequence of uFixedLevel
        # vector elements of the list are expected to be of non-uniform length, hence the requirement of collection
        # into a list, as opposed to a matrix or data frame 
        # it is assumed that FEX has already been exported
        #parLapply(cl, ulev, function(lev) which(FEX==lev))
      #}
      #environment(plyXTX) <- .GlobalEnv

      # construct indicator indices (pointers to fixed effect entries with value of 1)
      # construct high dimension effect indicators in parallel (the cost of memory export and
      # parallel overhead is offset by efficiencies only when dimension is sufficiently high; the
      # chosen value of 35 is an approximate boundary based on empirical observation, converting
      # compute durations of several minutes into approximately one minute)
      # for very high dimension (thousands of levels) inefficiencies of which() become apparent and the
      # alternative of ordering observation indices by level sorting levels demonstrates improved
      # performance (note that ordering is expensive when compared to parallelized which() operations,
      # but is compensated for by very fast match() operations)
      # empirical tests have shown parallel which() operations involving millions of levels require
      # many hours to complete, while order() and match() require approximately 1.5 minutes for 
      # 2.75 million levels in 25,000,000 observations (you read correctly)      
      #cl <- makePSOCKcluster(rep("localhost", nCoreXTX))
      kFE <- list()
      for(i in 1:length(uFixedLevel))
        # note the strict use of the sort.list(factor) method, which has been observed in all tests
        # to provide superior performance for all fixed effect dimensions
        if(nFixedLevel[i]>0) {
          # ordered levels with match() identify observation boundaries by level
          # do not index reference level observation levels
          #obsIndices <- intersect(order(fixedX[,i]), which(fixedX[,i]!=refLevel[i]))
          # sort.list, method="quick" is fast, but requires numerics, hence the conversion of
          # fixedX[,i] to a factor then to an integer
          # intersect eliminates reference level observation indices
          obsIndices <- intersect(sort.list(as.integer(factor(fixedX[,i])), method="quick", na.last=NA),
                                  which(fixedX[,i]!=refLevel[i]))
          # locate first observation index for each level
          firstObsIndex <- match(uFixedLevel[[i]], fixedX[obsIndices,i])
          # parse sets of observation indices for each non-reference level
          if(length(firstObsIndex)>1) {
            # first n-1 levels
            kFE[[i]] <- lapply(1:(length(firstObsIndex)-1),
                                  function(i) obsIndices[firstObsIndex[i]:(firstObsIndex[i+1]-1)])
            # last level
            kFE[[i]][[length(kFE[[i]])+1]] <- obsIndices[firstObsIndex[length(firstObsIndex)]:length(obsIndices)]
          } else {
            kFE[[i]] <- list(obsIndices)
          }
        #} else if(nFixedLevel[i]>35) {
          # parallel which() operations
          # creation of apparent duplicates of data frame vectors appears wasteful, but R does make
          # a separate copy of the referenced vectors unless modification is attempted, which we do not
        #  tcompute <- tcompute + proc.time()-t0
        #  t0 <- proc.time()
        #  FEX <- fixedX[,i]
          # export from local environment, in case operating within a function
        #  clusterExport(cl, "FEX", envir=environment())
        #  tmem <- tmem + proc.time()-t0
        #  t0 <- proc.time()
          # original parallel method - executed within a function causes export of entire local environment
          # we do not want to export anything beyond what has been explicitly exported above
          # kFE[[i]] <- parLapply(cl, uFixedLevel[[i]], function(lev) which(kFEX==lev))
          # lean method!
        #  kFE[[i]] <- plyXTX(cl, uFixedLevel[[i]])
        #} else {
        #  kFE[[i]] <- lapply(uFixedLevel[[i]], function(u) which(fixedX[,i]==u))
        }

      # release memory consumed by objects exported to parallel cores
      #stopCluster(cl)
      #rm(cl)
      gc()

      #############################################################################################################
      #### parse interaction specifications and standardize for analysis later
      #############################################################################################################

      # parse independent variables from interaction specification and identify as continuous (c) or fixed (f)
      # place continuous vars to left of relationship, if present
      # enumerate combinations of joint level values
      if(!is.null(interactionX)) {
        interactX <- do.call(rbind,
                       mapply(1:nrow(interactionX), SIMPLIFY=F,
                              FUN=function(i) {
                                    # classify both interacting variables as continuous or fixed effect
                                    # two possibilities for each
                                    vcl <- rbind(
                                             t(c("v"=interactionX[i,1],
                                                 "cl"=ifelse(interactionX[i,1] %in% colnames(contX), "c",
                                                      ifelse(interactionX[i,1] %in% colnames(fixedX), "f", NA)))),
                                             t(c("v"=interactionX[i,2],
                                                        "cl"=ifelse(interactionX[i,2] %in% colnames(contX), "c",
                                                             ifelse(interactionX[i,2] %in% colnames(fixedX), "f", NA)))))
                                    if(!is.na(vcl[1,2]) & !is.na(vcl[2,2])) {
                                      # swap v1 and v2 if v1 fixed and v2 continuous
                                      if(vcl[1,2]=="f" & vcl[2,2]=="c") {
                                        v <- vcl[1,]
                                        vcl[1,] <- vcl[2,]
                                        vcl[2,] <- v
                                      }
                                      # return col names, type, number of levels (product of number of levels of the cols)
                                      data.frame("x1"=vcl[1,1], "x2"=vcl[2,1], "class1"=vcl[1,2], "class2"=vcl[2,2],
                                                 "nCrossLevel"=ifelse(vcl[1,2]=="c", 1, nFixedLevel[which(colnames(fixedX)==vcl[1,1])]) *
                                                               ifelse(vcl[2,2]=="c", 1, nFixedLevel[which(colnames(fixedX)==vcl[2,1])]))
                                    } else {
                                      data.frame("x1"=NA, "x2"=NA, "class1"=NA, "class2"=NA, "nCrossLevel"=NA)
                                    }
                                  }))
      } else {
        interactX <- data.frame("x1"=character(), "x2"=character(), "class1"=character(), "class2"=character(), "nCrossLevel"=numeric())
      }

      nInteractX <- nrow(interactX)

      #########################################################################################################################
      # begin construction of X'X
      # when interactions not specified or no errors in specification
      #########################################################################################################################

      if(length(which(is.na(interactX[,"x1"])))==0) {

        #######################################################################################################################
        #### enumerate variables in design matrix (including expanded fixed effects), generate empty X'X
        #### construct column names (for continuous and fixed effects vars, interactions done later)
        #######################################################################################################################
    
        # define pointers to beginning and ending rows and columns in X'X for each variable class
        # constant
        jb0 <- 1
        # continuous vars (note that if continuous vars not specified, ending column is 1)
        jcontX <- c(2, 1+ncontX)
        # fixed effects
        # cycle through FE specs and offset column positions from continuous prior FE levels
        jfixedX <- t(mapply(1:length(nFixedLevel),
                       FUN=function(i) {
                             # compute final position of current FE (final continuous pos + prior FE levels +
                             # current FE levels)
                             pos1 <- jcontX[2] + sum(nFixedLevel[1:i])
                             # return beginning pos (final pos - number of levels + 1), ending pos
                             c(pos1 - nFixedLevel[i] + 1, pos1)
                           }))
        # interactions
        # cycle through interaction specs and offset column positions from fixed effects and prior interaction levels
        if(nInteractX>0)
          jinteractX <- t(mapply(1:nInteractX,
                            FUN=function(i) {
                                  # compute final position of current interaction (final fixed effect pos + prior interaction
                                  # levels + current interaction levels)
                                  pos1 <- jfixedX[nrow(jfixedX),2] + sum(interactX[1:i,"nCrossLevel"])
                                  # return beginning pos (final pos - number of levels + 1), ending pos
                                  c(pos1 - interactX[i,"nCrossLevel"] + 1, pos1)
                                }))

        # record X'X column dimension
        dimXTX <- 1+ncontX+sum(nFixedLevel)+sum(interactX[,"nCrossLevel"])

        # construct upper triangle of X'X
        # constant, continuous vars, and fixed effects

        # factors with greatest number of levels should appear last so that parallel ops (cores) are
        # distributed among longest lists
        XTX <- matrix(0, nrow=dimXTX, ncol=dimXTX, dimnames=list(1:dimXTX, 1:dimXTX))

        # col names (constant, continuous, and fixed effects)
        if(ncontX>0) {
          colnames(XTX)[1:(1+ncontX+sum(nFixedLevel))] <-
            c("b0", paste("x-", colnames(contX), sep=""),
              unlist(lapply(1:length(uFixedLevel), FUN=function(i) paste(colnames(fixedX)[i], "-", uFixedLevel[[i]], sep=""))))
        } else {
          colnames(XTX)[1:(1+sum(nFixedLevel))] <-
            c("b0",
              unlist(lapply(1:length(uFixedLevel), FUN=function(i) paste(colnames(fixedX)[i], "-", uFixedLevel[[i]], sep=""))))
        }

        #######################################################################################################################
        #### begin row-wise completion of upper triangle of X'X
        #### rows and cols corresponding to continuous and fixed effects variable products only
        #### cols corresponding to products of continuous/FE with interactions done in interaction section later
        #### use efficient (sparse) indicator columns for high dimension fixed effects and interactions
        #######################################################################################################################

        # row 1:  sums of constant, continuous, and fixed effects columns
        # b0,  sum of continuous vars, and count of fixed effect observations with non-0 level
        if(ncontX>0) {
          XTX[1,1:(1+ncontX+sum(nFixedLevel))] <- c( nX, apply(contX, 2, sum), unlist(lapply(kFE, function(k) lapply(k, length))))
        } else {
          XTX[1,1:(1+sum(nFixedLevel))] <- c( nX, unlist(lapply(kFE, function(k) lapply(k, length))))
        }

        # cross products of continuous vars with continuous variables and fixed effects
        # fixed effects indicators permit sum of var rows using pointers to indicator=1 rows
        if(ncontX>0)
          for(i in 1:ncontX) {
            i0 <- jcontX[1]+i-1
            # products of continuous vars with continuous vars
            # note that the cross product of individual columns (transpose of one times the other) within mapply
            # is significantly more efficient than multiplying the same colums as a set
            XTX[i0,i0:jcontX[2]] <- mapply(i:ncontX, FUN=function(j) t(contX[,i])%*%contX[,j])
            # products of continuous vars with fixed effects
            for(j in 1:length(kFE))
              XTX[i0,jfixedX[j,1]:jfixedX[j,2]] <- mapply(1:nFixedLevel[j], FUN=function(k) sum(contX[kFE[[j]][[k]],i]))
          }

        # cross products of fixed effects with fixed effects

        # define function to execute parallel operations, enabling it to be assigned to .GlobalEnv
        # see notes included in previous instance of function
        plyXTX <- function(cl, ife, ilev, jfe, njlev) {
          # note the following parApply parameters (example is for cross product of race and FY columns:
          # cl ....................... the cluster of cores
          # as.matrix(1:nFY) ......... a vector (but must be in matrix form) of FY indices
          #                            these are passed sequentially to the function as its first parameter
          # 1 ........................ pass FY indices row wise (sequentially down the single column we are supplying)
          # function(k, i) ........... function to execute for each k, i supplied
          # i ........................ the current loop index (current race)
          # note that a vector of length nFY is returned by parApply, one element for each k passed to the function
          # the jth element of this vector corresponds to the count of intersecting FY and jth bureau observations
          # that is, where both of current race and jth FY are 1
          # count of positions (observation indices) where both effects are 1
          # execute by fixed effect to take advantage of orthogonal properties of within effect columns
          parApply(cl, as.matrix(1:njlev), 1,
            function(jlev, ife, ilev, jfe) length(intersect(kFE[[ife]][[ilev]], kFE[[jfe]][[jlev]])), ife, ilev, jfe)
        }
        environment(plyXTX) <- .GlobalEnv

        # create parallel cluster and export indicator indices referenced in parallel routines
        # a little expensive, in time and memory, but compensated for through parallelization efficiencies
        tcompute <- tcompute + proc.time()-t0
        t0 <- proc.time()
        cl <- makePSOCKcluster(rep("localhost", nCoreXTX))
        # export from local environment, in case executing within a function
        clusterExport(cl, "kFE", envir=environment())
        tmem <- tmem + proc.time()-t0
        t0 <- proc.time()
        # construct each row, filling columns to the right
        for(ife in 1:length(kFE))
          for(ilev in 1:nFixedLevel[ife]) {
            i0 <- jfixedX[ife,1]+ilev-1
            # within effect columns are mutually exclusive
            # diagonal positions (i=j) are frequencies for corresponding level, all other cross products (i<>j) are 0
            XTX[i0,i0] <- length(kFE[[ife]][[ilev]])
            # cross product with all levels of remaining fixed effects
            jfe <- ife+1
            while(jfe<=length(kFE)) {
              XTX[i0,jfixedX[jfe,1]:jfixedX[jfe,2]] <-
                # local environment export version
                # parApply(cl, as.matrix(1:nFixedLevel[[jfe]]), 1,
                #   function(k, ife, ilev, jfe) length(intersect(kFE[[ife]][[ilev]], kFE[[jfe]][[k]])), ife, ilev, jfe)
                # replaced with non-export .GlobalEnv function
                plyXTX(cl, ife, ilev, jfe, nFixedLevel[[jfe]])
              jfe <- jfe+1
            }
          }

        gc()

        #######################################################################################################################
        #### continue with row-wise completion of upper triangle of X'X
        #### columns corresponding to continuous/FE row, interaction column cross products
        #### rows corresponding to interactions
        #### use efficient (sparse) indicator columns for high dimension fixed effects and interactions
        #######################################################################################################################

        if(nInteractX>0) {

          #####################################################################################################################
          #### construct names
          #####################################################################################################################

          colnames(XTX)[(2+ncontX+sum(nFixedLevel)):dimXTX] <-
            unlist(lapply(1:nInteractX,
              FUN=function(i)
                # if either of an interaction pair are continuous then v1 is continuous
                if(interactX[i,"class1"]=="c")
                  if(interactX[i,"class2"]=="c")
                    paste(interactX[i,"x1"], "-", interactX[i,"x2"], sep="")
                  else {
                    # point to fixed effect levels for second x
                    k <- which(colnames(fixedX)==interactX[i,"x2"])
                    paste(interactX[i,"x1"], "-", interactX[i,"x2"], "(", uFixedLevel[[k]], ")", sep="")
                  }
                else {
                  # both must be fixed
                  # point to fixed effect levels for both vars
                  kf1 <- which(colnames(fixedX)==interactX[i,"x1"])
                  kf2 <- which(colnames(fixedX)==interactX[i,"x2"])
                  paste(sort(rep(paste(interactX[i,"x1"], "(", uFixedLevel[[kf1]], ")-", sep=""), length(uFixedLevel[[kf2]]))),
                  paste(interactX[i,"x2"], "(", uFixedLevel[[kf2]], ")", sep=""), sep="")
                }))

          #################################################################################################################
          #### row 1: sums of element-wise products of interaction cols
          #################################################################################################################

          # type of sum depends on whether both vars are continuous, one continuous and one fixed, or both fixed
          for(i in 1:nInteractX) {
            v11 <- interactX[i,"x1"]
            v12 <- interactX[i,"x2"]
            if(interactX[i,"class1"]=="c")
              if(interactX[i,"class2"]=="c")
                # two continuous vars, sum element-wise products
                XTX[1,jinteractX[i,1]:jinteractX[i,2]] <- sum(contX[,v11]*contX[,v12])
              else {
                # continuous with a fixed effect, sum continuous elements corresponding to non-zero elements of fixed eff
                # identify index of fixed effect
                kf1 <- which(colnames(fixedX)==v12)
                XTX[1,jinteractX[i,1]:jinteractX[i,2]] <-
                  mapply(1:nFixedLevel[kf1], FUN=function(klev1) sum(contX[kFE[[kf1]][[klev1]],v11]))
              }
            else {
              # both must be fixed, count observations with 1 in both effects (intersection)
              # point to fixed effect levels for both vars
              kf1 <- which(colnames(fixedX)==v11)
              kf2 <- which(colnames(fixedX)==v12)
              # note that sums are produced for each FE 2 level nested within each FE 1 level
              XTX[1,jinteractX[i,1]:jinteractX[i,2]] <-
                unlist(lapply(1:nFixedLevel[kf1],
                  FUN=function(klev1) unlist(lapply(1:nFixedLevel[kf2],
                    FUN=function(klev2, klev1) length(intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])), klev1))))
            }
          }

          gc()

          #################################################################################################################
          #### cross products of continuous vars and interactions
          #### rows corresponding to continuous vars, cols corresponding to interactions
          #### these cols were left incomplete in continuous cross product section, above
          #################################################################################################################

          if(ncontX>0)
            for(i in 1:ncontX) {
              i0 <- jcontX[1]+i-1
              # compute cross-product of current continuous var with each interaction 
              # type of multiply operation depends on whether both interaction vars are continuous,
              # one continuous and one fixed, or both fixed
              for(inter1 in 1:nInteractX) {
                # identify interaction variables
                v21 <- interactX[inter1,"x1"]
                v22 <- interactX[inter1,"x2"]
                if(interactX[inter1,"class1"]=="c")
                  if(interactX[inter1,"class2"]=="c")
                    # interaction involves two continuous vars, sum element-wise products
                    XTX[i0,jinteractX[inter1,1]:jinteractX[inter1,2]] <-
                      sum(contX[,i]*contX[,v21]*contX[,v22])
                  else {
                    # interaction involves a continuous and a fixed effect
                    # sum product of current continuous var with interaction continuous var where fixed effect level is non-zero
                    # identify index of fixed effect (to filter continuous rows)
                    kf1 <- which(colnames(fixedX)==v22)
                    # PERFORMANCE OPPORTUNITY USING PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                    XTX[i0,jinteractX[inter1,1]:jinteractX[inter1,2]] <-
                      mapply(1:nFixedLevel[kf1],
                        FUN=function(klev1) sum(contX[kFE[[kf1]][[klev1]],i] * contX[kFE[[kf1]][[klev1]],v21]))
                  }
                else {
                  # both interaction vars fixed, sum continuous var with 1 in both fixed effects levels (intersection)
                  # point to fixed effect levels for both fixed vars
                  kf1 <- which(colnames(fixedX)==v21)
                  kf2 <- which(colnames(fixedX)==v22)
                  # note that sums are produced for each level of fixed effect 1 intersected with fixed effect 2
                  # PERFORMANCE OPPORTUNITY USING PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                  XTX[i0,jinteractX[inter1,1]:jinteractX[inter1,2]] <-
                    unlist(lapply(1:nFixedLevel[kf1],
                      FUN=function(klev1) unlist(lapply(1:nFixedLevel[kf2],
                        FUN=function(klev2, klev1)
                          sum(contX[intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]]),i]), klev1))))
                }
              }
            }

            gc()

          #################################################################################################################
          #### cross products of fixed effects and interactions
          #### rows corresponding to fixed effects, cols corresponding to interactions
          #### these cols were left incomplete in fixed effects cross product section, above
          #################################################################################################################

          # construct each row (level) of each fixed effect in sequence
          # note that execution reaches this point only when fixed effects have been declared
          # there is no need to test for length(kFE)>0)
          for(ife in 1:length(kFE))
            for(ilev in 1:nFixedLevel[ife]) {
              # XTX row pointer for current level of current fixed effect (row under construction)
              i0 <- jfixedX[ife,1]+ilev-1
              # observation indicators for current level of current fixed effect
              klev0 <- kFE[[ife]][[ilev]]
              # compute cross-product of current fixed effect var with each interaction 
              # type of multiply operation depends on whether both interaction vars are continuous,
              # one continuous and one fixed, or both fixed
              for(inter in 1:nInteractX) {
                # identify interaction variables
                v21 <- interactX[inter,"x1"]
                v22 <- interactX[inter,"x2"]
                if(interactX[inter,"class1"]=="c")
                  if(interactX[inter,"class2"]=="c")
                    # interaction involves two continuous vars, sum element-wise products where
                    # fixed level is non-zero
                    XTX[i0,jinteractX[inter,1]:jinteractX[inter,2]] <-
                      sum(contX[klev0,v21]*contX[klev0,v22])
                  else {
                    # interaction involves a continuous and a fixed effect
                    # sum continuous var elements where both fixed effect levels are non-zero
                    # identify index of interaction fixed effect (to intersect with current fixed effect)
                    kf1 <- which(colnames(fixedX)==interactX[inter,"x2"])
                    # PERFORMANCE OPPORTUNITY USING PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                    XTX[i0,jinteractX[inter,1]:jinteractX[inter,2]] <-
                      mapply(1:nFixedLevel[kf1],
                        FUN=function(klev) sum(contX[intersect(klev0, kFE[[kf1]][[klev]]),v21]))
                  }
                else {
                  # both interaction vars fixed, count elements in intersection of all three fixed effects
                  # point to fixed effect levels for both interaction fixed vars
                  kf1 <- which(colnames(fixedX)==interactX[inter,"x1"])
                  kf2 <- which(colnames(fixedX)==interactX[inter,"x2"])
                  # note that sums are produced for each level of fixed effect 1 intersected with fixed effect 2
                  # PERFORMANCE OPPORTUNITY USING PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                  XTX[i0,jinteractX[inter,1]:jinteractX[inter,2]] <-
                    unlist(lapply(1:nFixedLevel[kf1],
                      FUN=function(klev1) unlist(lapply(1:nFixedLevel[kf2],
                        FUN=function(klev2, klev1)
                          length(intersect(intersect(klev0, kFE[[kf1]][[klev1]]), kFE[[kf2]][[klev2]])), klev1))))
                }
              }
            }

            gc()

          #################################################################################################################
          #### cross products of interactions with interactions
          #### upper triangle rows and cols corresponding to cross products of interactions
          #### note that interaction cols of X'X were assigned in order of appearance in the interaction specification
          #### proceed in this order, computing cross products for an interaction with itself and trailing (to the right)
          #### interactions only
          #### note that interaction standardization, above, assures that interactions are of type
          #### continuous-continuous (cc), continuous-fixed (cf), or fixed-fixed (ff)
          #################################################################################################################

          # interaction cross products are the product of two interactions, each of which are the cross product of
          # two continuous vars, one continuous and one fixed effect, or two fixed effects vars  
          # this gives a total of nine possible interaction cross product types: (cc, cf, ff) X (cc, cf, ff), each with
          # a distinct optimal addressing method 

          # PERFORMANCE OPPORTUNITIES:
          # CONSOLIDATE COMMON SOLUTIONS (FOR INSTANCE, [c1Xc2]X[c3XF1] has solution structure of [c3Xf1]X[c1Xc2]
          # DO NOT EVALUATE OFF-DIAGONAL ELEMENTS (DIFFERING LEVELS) WHEN A FIXED EFFECT APPEARS MULTIPLE TIMES,
          # AS IN [c1Xf1]X[f1Xf2], SINCE WITHIN FIXED EFFECT VECTORS ARE ORTHOGONAL
          # PARALLELIZE!

          # inspect each declared interaction and execute optimal cross product methods involving current and
          # subsequent interactions
          # note that there exist nine permutations of interaction X interaction products (cc, cf, ff) X (cc, cf, ff)
          # identify mode of left and right operand then produce cross products exploiting intersections of non-zero
          # indices of all factors on either side of product
          # note that X'X rows to populate correspond to beginning and ending positions in jinteractX for left side
          # operand, while beginning and ending cols correspond to jinteractX values for right side operand

          # cycle through all interactions
          # note that, in order for processing to be at this point, the specified model includes interactions
          # therefore, there is no need to test nInteractX > 0 
          for(inter1 in 1:nInteractX)
            # cycle through current and all subsequent interactions
            # note that the first cycle involves the product of the left side interaction with itself
            for(inter2 in inter1:nInteractX) {

              # identify left (v11, v12) and right (v21 and v22) interacting vars
              v11 <- interactX[inter1,"x1"]
              v12 <- interactX[inter1,"x2"]
              v21 <- interactX[inter2,"x1"]
              v22 <- interactX[inter2,"x2"]

              # identify column boundaries of X'X elements to be updated
              # row boundaries are dependent on left side interaction type (cc, cf, ff)
              j0 <- jinteractX[inter2,1]
              j1 <- jinteractX[inter2,2]

              if(interactX[inter1,"class1"]=="c" & interactX[inter1,"class2"]=="c") {

                # interaction 1 (left operand) involves two continuous vars (cc)
                # use element-wise products of v11 and v12

                # identify X'X row offset
                i0 <- jinteractX[inter1,1]

                if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="c") {
                  # cc X cc
                  # interaction 2 (right operand) involves two continuous vars, use element-wise products of v21 and v22
                  XTX[i0,j0] <- sum(contX[,v11]*contX[,v12]*contX[,v21]*contX[,v22])
                } else if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="f") {
                  # cc X cf
                  # right side interaction contains a fixed effect
                  # use product of v11, v12, and v21 corresponding to non-zero indices of v22
                  # identify position in list of non-zero FE indices corresponding to v22
                  kf1 <- which(colnames(fixedX)==v22)
                  XTX[i0,j0:j1] <- mapply(1:nFixedLevel[[kf1]],
                                     FUN=function(klev1)
                                           # note that v11, v12, and v21 are continuous
                                           sum(contX[kFE[[kf1]][[klev1]],v11] *
                                                          contX[kFE[[kf1]][[klev1]],v12] *
                                                          contX[kFE[[kf1]][[klev1]],v21]))
                } else {
                  # cc X ff
                  # both vars in right side (v21, v22) interaction must be fixed effects since standardization places
                  # any continuous combination ahead of fXf
                  # use cross product of v11 and v12 in intersection of non-zero v21 and v22
                  # note that X'X interaction columns have been constructed in v1-level, v2-level sequence
                  # sums are produced for each combination of levels of v21 and v22
                  # PERFORMANCE OPPORTUNITY:  USE PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                  # point to fixed effect levels for both interaction fixed vars
                  kf1 <- which(colnames(fixedX)==v21)
                  kf2 <- which(colnames(fixedX)==v22)
                  XTX[i0,j0:j1] <-
                    unlist(lapply(1:nFixedLevel[kf1],
                             FUN=function(klev1)
                                   unlist(lapply(1:nFixedLevel[kf2],
                                            FUN=function(klev2, klev1) {
                                                  klev <- intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])
                                                  sum(contX[klev,v11]*contX[klev,v12])
                                                }, klev1))))
                }

              } else if(interactX[inter1,"class1"]=="c" & interactX[inter1,"class2"]=="f") {

                # interaction 1 involves one continuous var (v11) and one fixed var (v12) (cf)
                # use element-wise products involving v11 where v12 non-zero

                # identify non-zero FE indices corresponding to fixed var in left side interaction (v12)
                kf1 <- which(colnames(fixedX)==v12)

                # cycle through all levels of v12
                for(klev1 in 1:nFixedLevel[kf1]) {

                  # identify X'X row offset for current level of left side interaction FE (v12)
                  i0 <- jinteractX[inter1,1] + klev1 - 1

                  # update X'X cols based on right hand interaction style

                  if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="c") {
                    # cf X cc
                    # interaction 2 (right side) involves two continuous vars
                    # use element-wise products of v11, v21, and v22 where v12 non-zero
                    XTX[i0,j0] <- sum(contX[kFE[[kf1]][[klev1]],v11] *
                                                 contX[kFE[[kf1]][[klev1]],v21] *
                                                 contX[kFE[[kf1]][[klev1]],v22])
                  } else if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="f") {
                    # cf X cf
                    # right side interaction contains a fixed effect
                    # use cross product of v11 and v21 in intersection of non-zero v12 and v22
                    # PERFORMANCE OPPORTUNITIES:
                    #   USE PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                    #   DO NOT EVALUATE CROSS OFF DIAGONAL CROSS PRODUCTS WHEN LEFT AND RIGHT INTERACTIONS SHARE A FE
                    # identify position in list of non-zero FE indices for v22
                    kf2 <- which(colnames(fixedX)==v22)
                    XTX[i0,j0:j1] <- mapply(1:nFixedLevel[[kf2]],
                                       FUN=function(klev2) {
                                             # note that v11 and v21 are continuous
                                             klev <- intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])
                                             sum(contX[klev,v11] * contX[klev,v21])
                                           })

                  } else {
                    # cf X ff
                    # both vars in right side (v21, v22) interaction must be fixed effects since standardization places
                    # any continuous combination ahead of fXf
                    # use v11 elements in intersection of non-zero v12, v21, and v22
                    # note that X'X interaction columns have been constructed in v1-level, v2-level sequence
                    # sums are produced for each combination of levels of v21 and v22
                    # PERFORMANCE OPPORTUNITIES:
                    #   USE PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                    #   DO NOT EVALUATE CROSS OFF DIAGONAL CROSS PRODUCTS WHEN LEFT AND RIGHT INTERACTIONS SHARE A FE
                    # point to fixed effect levels for both right side interaction fixed vars
                    kf2 <- which(colnames(fixedX)==v21)
                    kf3 <- which(colnames(fixedX)==v22)
                    XTX[i0,j0:j1] <-
                      unlist(lapply(1:nFixedLevel[kf2],
                               FUN=function(klev2)
                                     unlist(lapply(1:nFixedLevel[kf3],
                                              FUN=function(klev3, klev2)
                                                    sum(contX[intersect(
                                                                intersect(kFE[[kf1]][[klev1]],
                                                                          kFE[[kf2]][[klev2]]),
                                                                          kFE[[kf3]][[klev3]]),v11])), klev2)))
                  }

                }

              } else if(interactX[inter1,"class1"]=="f" & interactX[inter1,"class2"]=="f") {

                # interaction 1 involves two fixed vars (v11 and v12) (ff)
                # use sums of intersections of all four vars

                # identify non-zero FE indices corresponding to fixed var in left side interaction (v12)
                kf1 <- which(colnames(fixedX)==v11)
                kf2 <- which(colnames(fixedX)==v12)

                # cycle through v11, v12 levels in X'X row order (all levels of v12 nested within all levels of v11)

                for(klev1 in 1:nFixedLevel[kf1])
                  for(klev2 in 1:nFixedLevel[kf2]) {

                    # identify X'X row offset for current level of v12 within current level of v11
                    i0 <- jinteractX[inter1,1] + (klev1-1)*nFixedLevel[kf2] + klev2 - 1

                    # identify joint non-zero elements of current levels of v11 and v12
                    klev0 <- intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])

                    # update X'X cols based on right hand interaction style

                    if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="c") {
                      # ff X cc
                      # interaction 2 (right side) involves two continuous vars
                      # use element-wise products of v21 and v22, indices from intersection of non-zero v11, v12
                      XTX[i0,j0] <- sum(contX[klev0,v21] * contX[klev0,v22])
                    } else if(interactX[inter2,"class1"]=="c" & interactX[inter2,"class2"]=="f") {
                      # ff X cf
                      # right side interaction contains a fixed effect
                      # use sum of v21 in intersection of non-zero v11, v12, and v22
                      # PERFORMANCE OPPORTUNITIES:
                      #   USE PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                      #   DO NOT EVALUATE CROSS OFF DIAGONAL CROSS PRODUCTS WHEN LEFT AND RIGHT INTERACTIONS SHARE A FE
                      # identify position in list of non-zero FE indices for v22
                      kf3 <- which(colnames(fixedX)==v22)
                      XTX[i0,j0:j1] <- mapply(1:nFixedLevel[[kf3]],
                                         FUN=function(klev3) sum(contX[intersect(klev0, kFE[[kf3]][[klev3]]),v21]))
                    } else {
                      # ff X ff
                      # both vars in right side (v21, v22) interaction are fixed effects
                      # use count of observations in inteaction of v11, v12, v21, and v22
                      # note that X'X interaction columns have been constructed in v1-level, v2-level sequence
                      # counts are produced for each combination of levels of v21 and v22
                      # PERFORMANCE OPPORTUNITIES:
                      #   USE PARLAPPLY IF INTERACTION LEVELS > OPTIMAL THRESHOLD
                      #   DO NOT EVALUATE CROSS OFF DIAGONAL CROSS PRODUCTS WHEN LEFT AND RIGHT INTERACTIONS SHARE A FE
                      # point to fixed effect levels for both right side interaction fixed vars
                      kf3 <- which(colnames(fixedX)==v21)
                      kf4 <- which(colnames(fixedX)==v22)
                      XTX[i0,j0:j1] <-
                        unlist(lapply(1:nFixedLevel[kf3],
                                 FUN=function(klev3)
                                       unlist(lapply(1:nFixedLevel[kf4],
                                                FUN=function(klev4, klev3)
                                                      length(intersect(intersect(klev0, kFE[[kf3]][[klev3]]),
                                                                       kFE[[kf4]][[klev4]])), klev3))))
                    }

                  }

              } # if inter1 cc, cf, ff

            } # for(inter2 in inter1:nInteractX)

        } # if(nInteractX>0)

        gc()

        ###################################################################################################################
        #### X'X upper triangle complete
        #### complete lower triangle (below diagonal) using transpose of upper
        ###################################################################################################################

        # note that lower triangle elements presently are all 0
        # note, also, the subtraction of the duplicate diagonal
        XTX <- XTX+t(XTX)-diag(x=diag(XTX), nrow=dimXTX)

        ###################################################################################################################
        #### construct X'Y
        ###################################################################################################################

        # intercept position
        XTY <- sum(Y)

        # transpose of continuous variables times Y
        # mapply of individual vector operations outperforms t(vectors)%*%Y
        if(ncontX>0)
          XTY <- c(XTY, mapply(1:ncontX, FUN=function(j) t(contX[,j])%*%Y))

        # fixed effects, sums of Y corresponding to FE=1 positions
        # compute t(FE indicator by level)%*%Y for each fixed effect
        # use previously constructed fixed effect indicator columns (kFE)
        # note that execution reaches this point only when fixed effects declared (no need to 
        # test for non-zero length(kFE)
        for(i in 1:length(kFE))
          XTY <- c(XTY, mapply(1:nFixedLevel[i], FUN=function(k) sum(Y[kFE[[i]][[k]]])))

        # interactions
        # compute t(interaction columns)%*%Y for each interaction
        # there exists three possible variable combinations: continuous X continuous (cc),
        # continuous X fixed (cf), and fixed X fixed (ff)
        # use previously constructed FE indicator columns (kFE) for cf and ff cases
        if(nInteractX>0)
          for(i in 1:nInteractX)
            if(interactX[i,"class1"]=="c" & interactX[i,"class2"]=="c") {
              # interaction of type cc, simple cross product
              XTY <- c(XTY, sum(contX[,interactX[i,"x1"]] * contX[,interactX[i,"x2"]] * Y))
            } else if(interactX[i,"class1"]=="c" & interactX[i,"class2"]=="f") {
              # interaction of type cf, return product of elements of continuous X and Y where
              # fixed effect is non-zero, evaluate for each level of fixed effect
              # identify position of fixed effect
              # note that, due to interaction standardization, FE always in var 2 of interaction for type cf
              kf1 <- which(colnames(fixedX)==interactX[i,"x2"])
              XTY <- c(XTY, mapply(1:nFixedLevel[kf1],
                              FUN=function(klev) sum(contX[kFE[[kf1]][[klev]],interactX[i,"x1"]] * Y[kFE[[kf1]][[klev]]])))
            } else {
              # interaction of type ff, return sum of Y elements where both fixed effects levels non-zero
              # note that one sum is produced for each level of FE2 nested within each level of FE1
              # identify positions of fixed effects
              kf1 <- which(colnames(fixedX)==interactX[i,"x1"])
              kf2 <- which(colnames(fixedX)==interactX[i,"x2"])
              XTY <- c(XTY, unlist(lapply(1:nFixedLevel[kf1], FUN=function(klev1)
                                     unlist(lapply(1:nFixedLevel[kf2], FUN=function(klev2, klev1)
                                              sum(Y[intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])]), klev1 )))) )
            }

        ###################################################################################################################
        #### coefficient solution
        ###################################################################################################################

        # solve X'X*Beta = X'Y, for Beta
        # simple, obsolete method:  beta <- solve(XTX, XTY, singular.ok=T)
        if(tolower(solMethod=="qr")) {
          # QR decomposition - will solve for parameter estimates even when XTX is singular
          # return QR solution along with rank
          # rank < XTX dimension => non-unique estimate solution (estimates for all but one
          # dependent column are returned as NA)
          qrXTX <- qr(XTX)
          rankXTX <- qrXTX$rank
          beta <- qr.coef(qrXTX, XTY)
        } else if(tolower(solMethod)=="chol") {
          # Cholesky decomposition using standard R chol() function
          # solve Ab=Y, where C'C=A
          # C'Cb=Y => C'z=Y (solve for z) and Cb=z (solve for b)
          # note that upper.tri applied before transpose
          chXTX <- chol(XTX)
          beta <- backsolve(chXTX, forwardsolve(chXTX, XTY, upper.tri=T, transpose=T), upper.tri=T)
          names(beta) <- colnames(XTX)
          rankXTX <- dimXTX  # need to adapt to rank of chol [attr(chXTX, "rank")]
        } else if(tolower(solMethod)=="chol-parallel") {
          # Cholesky decomposition using custom parallel Cholesky decomposition function
          # solve Ab=Y, where C'C=A
          # C'Cb=Y => C'z=Y (solve for z) and Cb=z (solve for b)
          # note that upper.tri applied before transpose
          chXTX <- choleskyDecomp(XTX)
          beta <- backsolve(chXTX, forwardsolve(chXTX, XTY, upper.tri=T, transpose=T), upper.tri=T)
          names(beta) <- colnames(XTX)
          rankXTX <- dimXTX  # need to adapt to rank of chol [attr(chXTX, "rank")]
        } else {
          status <- "Error:  Invalid solution method specified"
        }
        
        ###################################################################################################################
        #### estimate y values
        ###################################################################################################################

        # force collinear parameters (all but one) to 0
        betaLI <- beta
        betaLI[which(is.na(beta))] <- 0

        # constant
        yEst <- rep(betaLI[1], nX)

        # continuous variables
        if(ncontX>0)
          yEst <- yEst + as.matrix(contX)%*%betaLI[jcontX[1]:jcontX[2]]

        # fixed effects
        # FE existence already verified
        for(i in 1:length(kFE))
          for(j in 1:nFixedLevel[[i]])
            yEst[kFE[[i]][[j]]] <- yEst[kFE[[i]][[j]]] + betaLI[jfixedX[i,1]+j-1]

        # interactions
        # multiply interaction coefficients by computed interaction levels and add to Y
        # there exists three possible variable combinations: continuous X continuous (cc),
        # continuous X fixed (cf), and fixed X fixed (ff)
        # use previously constructed FE indicator columns (kFE) for cf and ff cases
        if(nInteractX>0)
          for(i in 1:nInteractX)
            if(interactX[i,"class1"]=="c" & interactX[i,"class2"]=="c") {
              # interaction of type cc, element-wise products
              yEst <- yEst + betaLI[jinteractX[i,1]] * contX[,interactX[i,"x1"]] * contX[,interactX[i,"x2"]]
            } else if(interactX[i,"class1"]=="c" & interactX[i,"class2"]=="f") {
              # interaction of type cf, update Y positions corresponding to non-zero FE elements
              # note that, due to interaction standardization, FE always in var 2 of interaction for type cf
              kf1 <- which(colnames(fixedX)==interactX[i,"x2"])
                for(klev1 in 1:nFixedLevel[kf1])
                  yEst[kFE[[kf1]][[klev1]]] <- yEst[kFE[[kf1]][[klev1]]] +
                                               betaLI[jinteractX[i,1]+klev1-1] *
                                               contX[kFE[[kf1]][[klev1]],interactX[i,"x1"]]
            } else {
              # interaction of type ff, update positions of Y corresponding to joint non-zero FE elements
              # cycle through each level of FE2 within each level of FE1
              # add FE1 X FE2 coefficient to joint non-zero observation positions
              # identify positions of fixed effects
              kf1 <- which(colnames(fixedX)==interactX[i,"x1"])
              kf2 <- which(colnames(fixedX)==interactX[i,"x2"])
                for(klev1 in 1:nFixedLevel[kf1])
                  for(klev2 in 1:nFixedLevel[kf2]) {
                    klev3 <- intersect(kFE[[kf1]][[klev1]], kFE[[kf2]][[klev2]])
                    yEst[klev3] <- yEst[klev3] + betaLI[jinteractX[i,1] + (klev1-1)*nFixedLevel[kf2] + klev2 - 1]
                  }
            }

        ###################################################################################################################
        #### estimate y-variance
        ###################################################################################################################

        yVar <- sum((Y-yEst)**2)/(nX-dimXTX)

        ###################################################################################################################
        #### X'X construction, beta solution, and y estimates are complete
        ###################################################################################################################

        status <- ""

      } else {

        status <- "ERROR (feXTX) - Invalid or unknown vectors in interaction specification"

      }


    } else {

      status <- "ERROR (feXTX) - No fixed effects specified in model"

    }

    # retain cluster if robust error option specified, sice fixed effect indicators used in parallel there
    if(tolower(estBetaVar)!="robust") {
      stopCluster(cl)
      rm(cl)
    }

    gc()

  }
  
  if(status=="") {

    ###################################################################################################################
    #### estimate parameter variances, if requested
    ###################################################################################################################

    if(tolower(estBetaVar)=="stdols") {

      # standard OLS estimates, assuming homoskedastic, iid errors
      # var(beta) = (variance of y estimates) * inverse(X'X)
      # verify nonsingularity of X'X, report NA otherwise
      if(rankXTX==dimXTX) {
        if(tolower(solMethod)=="qr") {
          vBeta <- yVar*diag(qr.solve(qrXTX))
        } else if(tolower(solMethod)=="chol") {
          # standard R chol2inv() function
          # note the use of prior cholesky decomposition
          vBeta <- yVar*diag(chol2inv(chXTX))
        } else if(tolower(solMethod)=="chol-parallel") {
          # diagonal of inverse using custom parallel Cholesky inverse function
          # note the use of prior cholesky decomposition
          vBeta <- yVar*cholInvDiag(chXTX)
        } else {
          vBeta <- rep(NA,dimXTX)
        }
      } else {
        vBeta <- rep(NA,dimXTX)
      }

    } else if(tolower(estBetaVar)=="robust" & nInteractX==0) {

      # robust, uncorrelated, heteroskedastic (independent, non-identically-distributed) errors
      # estimate robust variances (square of standard errors) as v(beta) = i(X'X)*X'uu'X*i(X'X)
      # where i() is the inverse function and u is the vector of observation response to predicted
      # value errors (y - yEstimate)
      # this derives from the expression for parameter variances from the normal OLS equaions
      # uu' forms an n X n diagonal error variance-covariance matrix with all off-diagonal elements
      # set to 0, which asserts the assumption of independent (uncorrelated) errors
      # but, since the diagonal elements are expected to be non-constant, this method models
      # parameter standard errors in a heteroskedastic (independent, but not identically distributed)
      # error setting 

      # since X'uu'X is symmetric, construct upper triangle only then copy transpose to lower triangle

      # verify nonsingularity of X'X, report NA otherwise
      if(rankXTX==dimXTX) {

        # create (0) X'uu'X matrix
        XuuX <- matrix(data=rep(0, dimXTX*dimXTX), nrow=dimXTX)

        # construct errors
        u <- (Y-yEst)**2

        # construct X'uu'X

        # row 1, constant, continuous vars, fixed effects
        if(ncontX>0) {
          XuuX[1,] <- c(sum(u), apply(as.matrix(1:ncontX), 1, function(i) sum(contX[,i]*u)),
                         unlist(lapply(1:length(kFE), function(ife)
                           unlist(lapply(1:nFixedLevel[ife], function(lev) sum(u[kFE[[ife]][[lev]]]))))))
        } else {
          XuuX[1,] <- c(sum(u), unlist(lapply(1:length(kFE), function(ife)
                                   unlist(lapply(1:nFixedLevel[ife], function(lev) sum(u[kFE[[ife]][[lev]]]))))))
        }

        # rows for continuous variables
        if(ncontX>0)
          for(i in 1:ncontX) {
            i0 <- jcontX[1]+i-1
            # products of continuous vars with continuous vars
            XuuX[i0,i0:jcontX[2]] <- mapply(i:ncontX, FUN=function(j) sum(contX[,i]*contX[,j]*u))
            # products of continuous vars with fixed effects
            for(ife in 1:length(kFE))
              XuuX[i0,jfixedX[ife,1]:jfixedX[ife,2]] <-
                apply(as.matrix(1:nFixedLevel[ife]), 1, 
                      function(lev) sum(contX[kFE[[ife]][[lev]],i]*u[kFE[[ife]][[lev]]]))
          }

        # rows for fixed effects
        # note that kFE has already been exported to parallel cluster (in XTX construction)
        clusterExport(cl, "u", envir=environment())
        # define function to execute parallel operations, enabling it to be assigned to .GlobalEnv
        # see notes included in previous instance of function
        plyVar <- function(cl, ife, ilev, jfe, njlev) {
          parApply(cl, as.matrix(1:njlev), 1,
            function(jlev, ife, ilev, jfe) sum(u[intersect(kFE[[ife]][[ilev]], kFE[[jfe]][[jlev]])]), ife, ilev, jfe)
        }
        environment(plyVar) <- .GlobalEnv
        # construct each row, filling columns to the right
        for(ife in 1:length(kFE))
          for(ilev in 1:nFixedLevel[ife]) {
            i0 <- jfixedX[ife,1]+ilev-1
            # within effect columns are mutually exclusive
            # diagonal positions (i=j) are 1 for observations in current effect  and level, others (i<>j) are 0
            XuuX[i0,i0] <- sum(u[kFE[[ife]][[ilev]]])
            # cross product with all levels of remaining fixed effects
            jfe <- ife+1
            while(jfe<=length(kFE)) {
              XuuX[i0,jfixedX[jfe,1]:jfixedX[jfe,2]] <-
                # local environment export version
                #apply(as.matrix(1:nFixedLevel[[jfe]]), 1,
                #  function(jlev, ife, ilev, jfe) sum(u[intersect(kFE[[ife]][[ilev]], kFE[[jfe]][[jlev]])]), ife, ilev, jfe)
                # parallel version
                plyVar(cl, ife, ilev, jfe, nFixedLevel[[jfe]])
              jfe <- jfe+1
            }
          }

        # complete lower triangle of XuuX
        # note that lower triangle elements presently are all 0
        # note, also, the subtraction of the duplicate diagonal
        XuuX <- XuuX+t(XuuX)-diag(x=diag(XuuX), nrow=nrow(XuuX))

        # finally, execute i(X'X)*X'uu'X*i(X'X)
        XTXi <- solve(XTX)
        vBeta <- diag(XTXi%*%XuuX%*%XTXi)

        gc()

      } else {

        vBeta <- rep(NA, dimXTX)

      }

    } else if(tolower(estBetaVar)=="cluster" & nInteractX==0) {

      # heteroskedastic, correlated within cluster, independent (uncorrelated) between cluster,
      # non-identically-distributed (heteroskedastic) errors
      # estimate clustered standard errors using v(beta) = i(X'X)*X'uu'X*i(X'X),
      # where i() implies inverse and uu' is the var-cov matrix of observation to estimate errors (y-yest)
      # this is the expression for var(beta) from ordinary least squares estimates of beta
      # on the assumption that X is non-stochastic and uu' is an estimate of error covariance and, further,
      # that covariance observed within clusters (groups) estimates that in the population while covariance between
      # groups is 0, uu' is altered such that all off-diagonal inter-group elements, that is uu'(ij) where
      # indices i and j belong to different groups, are set to 0
      # uu' has the form:
      #
      #     uu'(g1)    0       0    0 0 0 ...   0
      #        0    uu'(g2)    0    0 0 0 ...   0
      #        0       0    uu'(g3) 0 0 0 ...   0
      #                                   ...   0
      #        0       0       0      0   ... uu'(gk)
      #
      # where uu'(gi) is the n(gi) X n(gi) sub-matrix of the error covariace matrix corresponding to group i
      # note that with this construction of uu', X'uu'X is the sum (over all groups i) of (Xgi)'*uu'(gi)*Xgi,
      # where Xgi and uu'(gi) are the rows of X and sub-matrix of uu' corresponding to group i
      # [to see this, imagine multiplying X' by uu' and consider the effect of uu' elements of rows that do not
      # correspond to a particular group - each resulting column corresponds to products from a single group -
      # then multiplying that by X (each row belongs to a single group) gives a p X p sum of individual group products]

      # note that (Xg)'*ug = [(ug)'*Xg]', making (Xg)'*uu'(g)*Xg symmetric
      # also, (ug)'*Xg is a 1 X p vector, where p is the column dimension of X (the design matrix)
      # for each ID, construct v = (ug)'*Xg then use it to compute v'v = (Xg)'*uu'(g)*Xg

      # verify nonsingularity of X'X, report NA otherwise
      if(rankXTX==dimXTX) {

        # calculate errors
        u <- Y-yEst

        # generate fixed effect level indicators (pointers to levels for each observation)
        # note that reference levels contain NA
        xind <- mapply(1:ncol(fixedX), FUN=function(i) match(fixedX[,i], uFixedLevel[[i]]))
        colnames(xind) <- colnames(fixedX)

        # function to generate p X p sum of (Xg)'uu'(Xg) sub-matrices for pseudo IDs in specified range
        # X'uu'X = the sum of all individual (Xg)'uu'(Xg) sub-matrices
        # note that the upper triangle and diagonal of (symmetric) (Xg)'uu'(Xg) are returned
        XguuXg <- function(idboundk) {

          # create empty matrix to accumulate sum of group product matrices into
          XuuX <- matrix(data=rep(0, dimXTX*dimXTX), nrow=dimXTX)

          # Accumulate (Xg)'uu'(Xg) sub-matrices for individual pseudo IDs
          for(i in idbound[idboundk, 1]:idbound[idboundk, 2]) {

            # identify observations in current cluster ID group (all belonging to ith ID)
            # note that kg <- which(robustVarID==id) is very inefficient (hours with certain data sets) and is
            # replaced with the following
            # the first observation index for the current (ith) ID is in pidObsIndices[pidFirstObsIndex[i]]
            # the last observation index is in idObsIndices[idFirstObsIndex[i+1]-1] (this is the obs index
            # prior to the first index for ID i+1, the next ID)
            # the final index for the final ID is in the last position of idObsIndices since no ID follows
            # all observation indices between the first and final for an ID belong to that ID
            if(i<nid) {
              kg <- idObsIndices[idFirstObsIndex[i]:(idFirstObsIndex[i+1]-1)]
            } else {
              kg <- idObsIndices[idFirstObsIndex[i]:length(idObsIndices)]
            }

            # assemble (ug)'Xg, beginning with the intercept and continous terms
            # append 0s for indicator cols
            if(ncontX>0) {
              utXg <- c(sum(u[kg]), t(u[kg])%*%contX[kg,], rep(0, dimXTX-ncontX-1))
            } else {
              utXg <- c(sum(u[kg]), rep(0, dimXTX-1))
            }

            # get fixed effect level indices for each observation
            # note that reference levels have index NA
            kv <- matrix(xind[kg,], ncol=ncol(xind))

            # accumulate j=1..ng observation errors into positions corresponding to each FE level index
            # NA (reference level) indices are excluded
            # note that jfixedX[kv2,1]-1 are pointers into X for the first level of each fixed effect
            for(j in 1:length(kg)) {
              kv2 <- which(!is.na(kv[j,]))
              kv2 <- jfixedX[kv2,1]-1+kv[j,kv2]
              utXg[kv2] <- utXg[kv2] + u[kg][j]
            }

            # accumulate sub-matrix products
            # it is tempting to execute XuuX <- XuuX + utXg%*%t(utXg), but this is computationally expensive
            # (adds approximately 0.08 seconds per ID in trials) and involves products of which over 99% involve a factor of 0
            # it also calculates products in both upper and lower triangles of XuuX, when XuuX is known to be symmetric
            # alternative:  use products of non-zero entries of utXg and accumulate by index into appropriate rows and columns of XuuX
            # complete the upper triangle only, it will be copied to the lower triangle later
            # get non-zero indices
            kv <- as.vector(which(utXg!=0))
            # permutate all non-zero indices
            kv <- cbind(sort(rep(kv, length(kv))), kv)
            # retain non-zero indices (ij) where j>=i
            kv <- kv[kv[,2]>=kv[,1],]
            # accumulate non-zero (ug)'Xg index products into cells of XuuX corresponding to those indices (upper triangle)
            XuuX[kv] <- XuuX[kv] + utXg[kv[,1]]*utXg[kv[,2]]

          }

          return(XuuX)

        }

        # list IDs
        uid <- sort(unique(data[,robustVarID]))

        # construct boundaries of ID indices, one for each core
        nid <- length(uid)
        nidblock <- as.integer(nid/nCoreVar)
        idbound <- cbind(0:(nCoreVar-1)*nidblock+1, (1:nCoreVar)*nidblock) 
        idbound[nCoreVar,2] <- nid

        # generate observation pointers by ID
        # note that this, somewhat complicated method, replaces use of which(robustVarID==id) in the
        # XguuXg function due to its extreme inefficiency (adds hours to execution time with millions of observations)
        # the following method generates a matrix of observation indices for each ID in ID order
        # list observation indices in ID order
        # ordering of IDs requires approximately one minute (25,000,000 observations)
        # match requires approximately three seconds (fast, huh?)
        idObsIndices <- order(data[,robustVarID])
        # locate first observation index for each ID
        idFirstObsIndex <- match(uid, data[idObsIndices,robustVarID])

        # create parallel cluster
        # release existing cluser, if present
        tcompute <- tcompute + proc.time()-t0
        t0 <- proc.time()
        if(exists("cl"))
          stopCluster(cl)
        cl <- makePSOCKcluster(rep("localhost", nCoreVar))
        # export objects referenced in parallel instructions
        # from local environment, in case executing within a function
        clusterExport(cl, c("idbound", "xind", "dimXTX", "u", "ncontX", "jfixedX", "contX", "idObsIndices",
                            "idFirstObsIndex", "nid"), envir=environment())
        tmem <- tmem + proc.time()-t0
        t0 <- proc.time()

        gc()
        
        # define function to execute parallel operations
        # this enables assigning the global environment as its env, which avoids parLapply's
        # insistence on exporting the entire local environment (it never exports .GlobalEnv,
        # why it exports anything is a puzzle since environment objects cannot be referenced
        # in parallel functions anyway) 
        # this is important if executing within a function with a local environment
        # since it prevents export of potentially large objects
        plyVar <- function(cl, idbound, XguuXg) parLapply(cl, 1:nrow(idbound), XguuXg)
        environment(plyVar) <- .GlobalEnv

        # call, in parallel, the ID-sub-matrix function
        # parLapply returns a list, the ith element containing X'uu'X sub-matrix for ID (group) i
        # Reduce element-wise sums the sub-matrices, giving the complete (upper triangle of) X'uu'X
        XuuX <- Reduce("+", plyVar(cl, idbound, XguuXg))

        # release parallel cores and memory
        stopCluster(cl)
        rm(cl)
        gc()

        # complete lower triangle of XuuX
        # note that lower triangle elements presently are all 0
        # note, also, the subtraction of the duplicate diagonal
        XuuX <- XuuX+t(XuuX)-diag(x=diag(XuuX), nrow=nrow(XuuX))

        XTXi <- solve(XTX)
        vBeta <- diag(XTXi%*%XuuX%*%XTXi)

      } else {

        vBeta <- rep(NA, dimXTX)

      }
      
    } else if(tolower(estBetaVar)=="i(x'x)x'uu'xi(x'x)" & nInteractX==0) {

      # this option is added as a test of interior X and u=(y-yEst) structure
      # analytically, X'u = X'(y-yEst) = X'y - X'X*beta = X'[X*inverse(X'X)X'y) = X'y - X'y = 0
      # so the following should result in a vector of length p (one element for each beta) each
      # containing near zero entries
      # note that X'uu'X is (n X n) symmetric, resulting from the product of X'u (n X 1) and u'X (1 X n)
      # method:  construct X'u then assign product of elements i and j of X'u to cell ij of X'uu'X
      # since X'uu'X is symmetric, construct upper triangle only then copy transpose to lower triangle

      # verify nonsingularity of X'X, report NA otherwise
      if(rankXTX==dimXTX) {

        # create empty X'u vector and X'uu'X matrix
        Xu <- vector("numeric", dimXTX)
        XuuX <- matrix(data=rep(0, dimXTX*dimXTX), nrow=dimXTX)

        # construct errors
        u <- Y-yEst

        # construct X'u
        # position 1 of X'u, the sum of u elements since col 1 of X = 1
        Xu[1] <- sum(u)

        # continuous vars
        if(contX>0)
          Xu[jcontX[1]:jcontX[2]] <- t(contX)%*%u

        # fixed effect columns contain either 0 or 1 making cross products, for a given level,
        # the sum of u elements corresponding to non-zero elements of the fixed effect level
        # cycle through all fixed effects, cross multiplying each level
        for(i in 1:length(kFE))
          Xu[jfixedX[i,1]:jfixedX[i,2]] <- mapply(1:nFixedLevel[i], FUN=function(j) sum(u[kFE[[i]][[j]]]))

        # X'u is composed, now cross multiply X'u and u'X
        # compose ij pairs, where j>=i (upper triangle with diagonal)
        ij <- cbind(sort(rep(1:dimXTX, dimXTX)), 1:dimXTX)
        ij <- ij[which(ij[,2]>=ij[,1]),]

        # copy X'u cross products to X'uu'X
        XuuX[ij] <- Xu[ij[,1]]*Xu[ij[2]]

        # complete lower triangle of XuuX
        # note that lower triangle elements presently are all 0
        # note, also, the subtraction of the duplicate diagonal
        XuuX <- XuuX+t(XuuX)-diag(x=diag(XuuX), nrow=nrow(XuuX))

        # finally, execute i(X'X)*X'uu'X*i(X'X)
        XTXi <- solve(XTX)
        vBeta <- diag(XTXi%*%XuuX%*%XTXi)

        gc()

      } else {

        vBeta <- rep(NA, dimXTX)

      }
        
    } else {

      vBeta <- rep(NA, dimXTX)

    }

    names(vBeta) <- colnames(XTX)
    status <- ""

  }

  ###################################################################################################################
  #### processing complete
  #### clean up and report results
  ###################################################################################################################

  if(status=="") {
  
    if(exists("cl")) {
      stopCluster(cl)
      rm(cl)
    }
    gc()
    tcompute <- tcompute+proc.time()-t0
    #print(paste("compute time (user, sys, elapsed):  ", round(tcompute["user.self"], 2), ", ", round(tcompute["sys.self"], 2), ", ", round(tcompute["elapsed"], 2), sep=""))
    #print(paste("mem export time (user, sys, elapsed):  ", round(tmem["user.self"], 2), ", ", round(tmem["sys.self"], 2), ", ", round(tmem["elapsed"],  2), sep=""))
    list("status"=status, "Y"=fparams$Y, "contX"=fparams$contX, "fixedX"=fparams$fixedX, "refLevel"=refLevel,
         "interactionX"=interactionX, "nX"=nX, "beta"=beta, "vBeta"=vBeta, "estBetaVar"=estBetaVar,
         "XTX"=XTX, "dimXTX"=dimXTX, "rankXTX"=rankXTX, "XTY"=XTY, "yEst"=as.vector(yEst), "yVar"=yVar,
         "time"=rbind(tcompute, tmem))

  } else {
  
    gc()
    list("status"=status)
    
  }
  
}