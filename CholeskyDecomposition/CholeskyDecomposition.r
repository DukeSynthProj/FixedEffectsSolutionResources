options(max.print=1000)      # number of elements, not rows
options(stringsAsFactors=F)
options(scipen=999999)

library(Rcpp)

#######################################################################################
#### Cholesky decomposition
#### method utilizing array of pointers to array pointers for large matrix (beyond 
#### 10,000 X 10,000)
#######################################################################################

# Cholesky decomposition
cSource <- "
  #include <Rcpp.h>
  #include <omp.h>
  // [[Rcpp::plugins(openmp)]]
  // following prefixes Rcpp type with Rcpp::
  using namespace Rcpp;
  // create array of p vector pointers, each addressing one matrix row
  // this compensates for the following issues when the input matrix is large
  // local arrays (created within function) of size greater than [2047,2047] cause
  // crash on execution
  // static arrays of size greater than approximately [10000,10000] cause
  // compilation errors, assembler messages are returned:  value of x too large
  // for field of 4 bytes, possibly due to, say, 20,000 X 20,000 =
  // 400,000,000 * 8 (double) = 3.2 billion, which is greater than unsigned long
  // capacity (it is assumed that assembler instructions must be generated to
  // address matrix positions offset from position 0)
  // declaring a matrix with long long constant, as in static double
  // R1[20000LL][20000LL]; also produces compiler 'too large' message
  // it is assumed that the input matrix, X, is square, symmetric, and positive-
  // definite
  // [[Rcpp::export]]
  NumericMatrix choleskyDecomp(NumericMatrix X) {
    // use long indices for large addressing, but signed since i is decremented to -1
    long p=X.ncol();
    long i, j, k;
    double ssqRii, s;
    // create p-length arrays of pointers, each addressing one p-length matrix row
    double **R = new double *[p];
    // create array of matrix row, one pointer for each row
    for(i=0; i<p; i++) {
      R[i] = new double [p];
      for(j=0; j<p; j++)
        R[i][j]=0;
    }
    // note that, where the analytical matrix equations solve for Rij by accumulating
    // prior computed elements in order of i, the following algorithm places computed
    // values for analytical Rij in R[j][i], causing accumulation accumulation along
    // columns, taking advantage of the row-major matrix storage of C (successive
    // columns are in contiguous memory)
    // in empirical testing, this improves performance by a factor of three
    // note, also, the corresponding return of R', placing Rji in Rij
    // compute row 0, it is trivial and referenced in subsequent rows
    R[0][0]=sqrt(X(0,0));
    for(j=1; j<p; j++)
      R[j][0]=X(0,j)/R[0][0];
    // compute rows 1 through p-1
    for(i=1; i<p; i++) {
      // accumulate sum of squared Rij elements, needed for R(i+1,i+1)
      ssqRii=0;
      for(k=0; k<i; k++)
        ssqRii+=R[i][k]*R[i][k];
      // compute diagonal element of current row
      R[i][i]=sqrt(X(i,i)-ssqRii);
      // compute off diagonal elements, Rij for j>i
      // recall R is upper diagonal by definition
      if(i<p-1)
        #pragma omp parallel for private(j, k, s)
        for(j=i; j<p; j++) {
          // accumulate prior computed Rki*Rkj products
          s=0;
          for(k=0; k<i; k++)
            s+=R[i][k]*R[j][k];
          // compute Rij using Xij and sum of prior Rki*Rkj products
          R[j][i]=(X(i,j)-s)/R[i][i];
        }
    }
    // copy R to array to be output
    // this is necessary since R is constructed as an array of independent pointers
    // to rows of the decomposition
    NumericMatrix Rout(p,p);
    std::fill(Rout.begin(), Rout.end(), 0);
    for(i=0; i<p; i++)    
      for(j=i; j<p; j++)
        Rout(i,j)=R[j][i];
    // release memory allocated to dynamic rows
    for(i=0; i<p; i++)
      delete [] R[i];
    // release memory allocated to arrays of row pointers
    delete [] R;
    return(Rout);
  }"

# compile and create .dll
# cacheDir contains .dll and R instructions for loading and creating function from .dll
sourceCpp(code=cSource, rebuild=T, showOutput=T, cacheDir=getwd(), cleanupCacheDir=F)





#######################################################################################
#### diagonal of inverse of Cholesky decomposition
#### method utilizing array of pointers to array pointers for large matrix (beyond 
#### 10,000 X 10,000)
#### observed performance is identical to the contiguous R[p][p], iX[p][p]
#### implementation with X'X of p=1,228
#######################################################################################

# diag of Cholesky inverse
cSource <- "
  #include <Rcpp.h>
  #include <omp.h>
  // [[Rcpp::plugins(openmp)]]
  // following prefixes Rcpp type with Rcpp::
  using namespace Rcpp;
  // create array of p vector pointers, each addressing one matrix row
  // this compensates for the following issues when the input matrix is large
  // local arrays (created within function) of size greater than [2047,2047] cause
  // crash on execution
  // static arrays of size greater than approximately [10000,10000] cause
  // compilation errors, assembler messages are returned:  value of x too large
  // for field of 4 bytes, possibly due to, say, 20,000 X 20,000 =
  // 400,000,000 * 8 (double) = 3.2 billion, which is greater than unsigned long
  // capacity (it is assumed that assembler instructions must be generated to
  // address matrix positions offset from position 0)
  // declaring a matrix with long long constant, as in static double
  // R1[20000LL][20000LL]; also produces compiler 'too large' message
  // it is assumed that the input matrix, X, is square, symmetric, and positive-
  // definite
  // [[Rcpp::export]]
  NumericVector cholInvDiag(NumericMatrix R0) {
    // use long indices for large addressing, but signed since i is decremented to -1
    long p=R0.ncol();
    long i, j, k;
    double s;
    // create p-length arrays of pointers, each addressing one p-length matrix row
    // a copy of R0 is placed in the R1 pointers, one row per pointer
    // A is the computed inverse of the matrix corresponding to R0 (R0'R0=X)
    // the vector d will contain the parsed diagonal of A and is returned
    double **R1 = new double *[p];
    double **A = new double *[p];
    NumericVector d(p);
    // create matrix row arrays, one for each row pointer
    // divide rows of (R) by diagonals (sets diag to 1)
    // eliminating need of division by R[i,i] in
    // constructing elements in row i+1
    for(i=0; i<p; i++) {
      R1[i] = new double [p];
      A[i] = new double [p];
      for(j=0; j<p; j++)
        R1[i][j]=R0(i,j)/R0(i,i);
      // diagonals of right hand matrix are squared
      // reciprocals of original R diagonals
      // they are the initial values of the inverse
      A[i][i]=1/R0(i,i)/R0(i,i);
    }
    // construct upper rows of the inverse in reverse
    // order, copy to transpose positions for use by
    // following row construction
    // note that rows cannot be constructed in parallel
    // since prior rows are necessary to construct
    // subsequent rows
    //omp_set_num_threads();
    for(i=p-2; i>=0; i--) {
      // construct columns in parallel since elements of
      // rows are referenced
      #pragma omp parallel for private(j, k, s)
      for(j=i+1; j<p; j++) {
        s=0;
        // following iterative loop is the computationally
        // intensive segment of this algorithm
        // note that, although the analytical matrix
        // equations multiply R[i,k]' by A[k,j], it
        // is observed that A is symmetric which allows
        // the equivalent operation R[i,k]' * A[j,k],
        // which is considerably more efficient as
        // implemented by the C compiler, since C stores
        // a matrix in row dominant order, making elements
        // adjacent in memory with respect to k (note that
        // each increment of j corresponds to an in memory
        // distance of p*8, where 8 is the number of bytes
        // per double floating point element)
        for(k=i+1; k<p; k++)
          s+=R1[i][k]*A[j][k];
        A[i][j]=-s;
        // copy to transpose position
        A[j][i]=-s;
      }
      // compute diagonal element on current row
      for(j=i+1; j<p; j++)
        A[i][i]=A[i][i]-R1[i][j]*A[i][j];
    }
    // copy all diagonal elements of computed inverse to result vector
    for(i=0; i<p; i++)
      d(i)=A[i][i];
    // release memory allocated to dynamic rows
    for(i=0; i<p; i++) {
      delete [] R1[i];
      delete [] A[i];
    }
    // release memory allocated to arrays of row pointers
    delete [] R1;
    delete [] A;
    return(d);
  }"

# compile and create .dll
# cacheDir contains .dll and R instructions for loading and creating function from .dll
sourceCpp(code=cSource, rebuild=T, showOutput=T, cacheDir=getwd(), cleanupCacheDir=F)








