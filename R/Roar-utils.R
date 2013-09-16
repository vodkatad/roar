getFisher <- function(counts) {
   # right_pre right_post left_pre left_post
   mat <- matrix(counts, ncol=2, byrow=FALSE)
   # Reminders:
   #colnames(mat) <- c('right', 'left')
   #rownames(mat) <- c('PRE', 'POST')
   f <- fisher.test(as.table(mat), alternative="two.sided")
   return(f$p.value)
}

checkStep <- function(rds, neededStep) {
   if (rds@step > neededStep) { # Already done, not repeating. Give a reset method?
      warning("Not repeating a step already done")
      return(c(FALSE, rds)) # We hope for promotion.
   } else if (rds@step < neededStep) { # Something is missing.
      if (neededStep != rds@step) {
         if (neededStep == 1) {
            warning("Automatically calling countPrePost")
            rds <- countPrePost(rds, FALSE)    
         } else if (neededStep == 2) {
            warning("Automatically calling computeRoars")
            rds <- computeRoars(rds)
         } else if (neededStep == 3) {
            warning("Automatically calling computePvals")
            rds <- computePvals(rds)
         }
      }
   } 
   return(c(TRUE, rds))
}

meanAcrossAssays <- function(assays, wantedColumns) {
   # Is the conversion to dataframe slow? Is there a more efficient way without
   # changing the SE/Assays structure?
   wantedCols <- lapply(assays, function(x) { x[,wantedColumns] } )
   return(rowMeans(as.data.frame(wantedCols)))
}