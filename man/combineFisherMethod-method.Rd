\name{combineFisherMethod}
     \alias{combineFisherMethod}
     \title{Private/inner/helper method to perform Fisher test}
     \description{
       This method \bold{should not} be used by package users. Given a numerical vector
       of length 4 it will perform a Fisher test and return the p-value for the two.sided
       test. Non-integer values will be rounded.
     }
     \usage{
      combineFisherMethod(pvals)
     }
     \arguments{
       \item{pvals}{
           A numerical vector of length 4.
        }
      }
      \value{
         The pvalue for the two.sided Fisher test.
      }
\keyword{combineFisherMethod}