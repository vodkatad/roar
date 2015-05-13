getFisher <- function(counts) {
   # treatment_pre treatment_post control_pre control_post
   mat <- matrix(counts, ncol=2, byrow=FALSE)
   # Reminders:
   #colnames(mat) <- c('treatment', 'control')
   #rownames(mat) <- c('PRE', 'POST')
   f <- fisher.test(as.table(mat), alternative="two.sided")
   return(f$p.value)
}

combineFisherMethod <- function(pvals) {
   # Computes the combined pvalue from a list of independent tests resulting pvalues
   # following the fisher method.
   # http://en.wikipedia.org/wiki/Fisher%27s_method
   chisq <- -2*sum(log(pvals))
   return(1-pchisq(chisq, df = 2*length(pvals)))
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

# XXX we need to have the seqlengths to have a correct last intron. TODO
# maybe we could avoid doing this...we need to have faith in the gtf tough
getApaGenesFractionsPlusStrand <- function(geneGr, apaGr)
{
   introns <- gaps(geneGr)
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how. 
   # New idea: we combine exons and introns adding them mcols to their kind
   # and with info about being overlappant or not.
   # Then we run over this GRanges that should have all the needed info
   # to generate a GRanges res with intervals and maybe also a list of 
   # pre/post indexes for every APA?
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'
   mcols(introns)$type <- 'i'   
   whole <- sort(c(geneGr, introns))
   hits <- findOverlaps(whole, apaGr)
   if (!all(countSubjectHits(hits) == 1)) {
      stop("Error: a given apa does not overlap its gene")
   }
   foundov <- length(unique(whole[queryHits(hits)]))
   mcols(whole)$overlap <- rep(FALSE, length(whole))   
   mcols(whole[unique(queryHits(hits))])$overlap <- rep(TRUE, foundov)
   if (whole[1]$type != 'i' || whole[1]$overlap) {
      stop("Error in the given gene structure or overlapping apas")
   }
   # Growing objects sin.
   begins <- ()
   ends <- ()
   begin <- -1
   # Write logic and invariant here TODO.
   for i in 2:length(whole) {
      if (whole[i]$type=='e') {
         if (whole[i]$overlap || (i+1<=length(whole) && whole[i+1]$overlap)) {
            if (begin != -1) {
               begins <- c(begins, begin)
               ends <- c(ends, end(whole[i-1]))
            }
            begin <- start(whole[i])
         }
      }
      if (whole[i]$overlap) {
         overlapping_apas <- apasGr[subjectHits(hits[queryHits(hits)==i])]
         for (j in 1:overlapping_apas) {
            begins <- c(begins, begin)
            ends <- c(ends, start(overlapping_apas[j]))
            begin <- end(overlapping_apas[j])
         }
      }
   }
   # We could be ok like this (last APA inside an exon, we stop there)
   # or we could have a last APA after the last exon that we want to consider.
   if (an APA outside hits) {
      begins <- c(begins, begin)
      ends <- c(ends, start(the last apa))
   }
}