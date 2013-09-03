\name{cores}
     \alias{cores}
     \title{Method to check how many cores are used by a roar analysis - right now not useful}
     \description{
         Right now always returns 1 as long as multi-core support has to be implemented.
     }
     \usage{
      cores(rds)
     }
     \arguments{
       \item{rds}{
           A \code{\linkS4class{RoarDataset}} object.}
      }
      \value{
         The number of cores used by this roar analisys.
      }
\keyword{checkStep}