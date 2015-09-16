\name{meanAcrossAssays}
     \alias{meanAcrossAssays}
     \title{Private/inner/helper method to get average counts across samples}
     \description{
       This method \bold{should not} be used by package users. It gets average counts
       for "pre" or "post" portions (depending on the wantedColumns argument)
       given the list of assays for one of the two conditions.
     }
     \usage{
      meanAcrossAssays(assays, wantedColumns)
     }
     \arguments{
       \item{assays}{
           A list of matrices/dataframes.
        }
        \item{wantedColumns}{
           The name of the columns ("pre" or "post") whose means should be computed.
           Average will be calculated on the corresponding rows of the list of matrices/dataframe,
           working on the given column.
        }
      }
      \value{
         The pvalue for the two.sided Fisher test.
      }

\keyword{meanAcrossAssays}