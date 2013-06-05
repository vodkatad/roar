get_fisher <- function(counts) {
   # right_pre right_post left_pre left_post
   mat <- matrix(counts, ncol=2)
   # Reminders:
   #colnames(mat) <- c('right', 'left')
   #rownames(mat) <- c('PRE', 'POST')
   f <- fisher.test(as.table(mat), alternative="two.sided")
   return(f$p.value)
}

checkStep <- function(rds, neededStep) {
   if (rds@step > neededStep) { # Already done, not repeating. Give a reset method?
      # Warning
      return(c(FALSE, rds)) # We hope in promotion.
   } else if (rds@step < neededStep) { # Something is missing
      while (neededStep != rds@step) {
         switch(neededStep,
            1= { rds <- countPrePost(rds) },
            2= { rds <- computeRoars(rds) },
            3= { rds <- computePvals(rds) }
         )
      }
   } 
   return(c(TRUE, rds))
}