\name{combineFisherMethod}
     \alias{combineFisherMethod}
     \title{Private/inner/helper method to combine pvalues of independent test with the Fisher method}
     \description{
       This method \bold{should not} be used by package users. Given a numerical vector
       of pvalues, which should be obtained from independent tests on the same null hypothesis, 
         this will give the combined pvalue following the Fisher method.
     }
     \usage{
      combineFisherMethod(pvals)
     }
     \arguments{
       \item{pvals}{
           A numerical vector with pvalues of independent tests on the same H0.
      }}
      \value{
         The combined pvalue given by the Fisher method.
      }
\keyword{combineFisherMethod}