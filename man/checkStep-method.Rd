\name{checkStep}
     \alias{checkStep}
     \title{Private/inner/helper method to check the order of the invoked analysis methods}
     \description{
       This method \bold{should not} be used by package users. It gets an rds object and a
       required number of analysis step and, if possible, calls the requested method to reach
       that step. It returns the object and a logical value that tells if the analysis can go on.
     }
     \usage{
      checkStep(rds, neededStep)
     }
     \arguments{
       \item{rds}{
           A \code{\linkS4class{RoarDataset}} object.}
       \item{neededStep}{
           The analysis step where rds should be/arrive.    
       }
      }
      \value{
         A list containing a logical that shows if the needed step could be reached
         with rds and the object at the requested step. Check step won't repeat
         a step already done and the logical value will be FALSE in this case (and rds
         won't be returned modified).
      }
\keyword{checkStep}