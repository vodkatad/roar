get_fisher <- function(counts) {
   # right_pre right_post left_pre left_post
   mat <- matrix(counts, ncol=2)
   # Reminders:
   #colnames(mat) <- c('right', 'left')
   #rownames(mat) <- c('PRE', 'POST')
   f <- fisher.test(as.table(mat), alternative="two.sided")
   return(f$p.value)
}