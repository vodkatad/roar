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

getApaGenesFractions <- function(geneGr, apaGr)
{
   introns <- gaps(geneGr)
   #mcols(geneGr) <- NULL
   #all_int <- c(geneGr, introns)
   #info <- c()
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how. subsetByOverlaps(gr,grl) could be an idea.
   #apa_n <- 1
   #for i in 1:length(all_int) {
   #   if
   #}
   hits_ex <- findOverlaps(geneGr, apaGr)
   hits_int <- findOverlaps(intronsGr, apaGr)
   last_int <- length(introns) # APAs falling here will be outsiders
   if (!all(countSubjectHits(hits_ex)+countSubjectHits(hits_int) == 1)) {
      stop "Error: a given apa does not overlap its gene"
   }
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how.
   chr <- "chr" # get from them always it
   strand <- "strand"
   GRanges res <- GRanges()
   begin <- start(head(geneGr, n=1))
   # Type is the kind of beginning of this interval. For
   # each APA we define two intervals (but not for outsiders).
   # We are in the circle of growing objects. FIXME.
   # Reason more about which intervals are needed. We need to
   # have beginning of exons to define pre/post and not only APA defined
   # portions as we are doing here.
   for i in 1:length(apaGr) {
      if (length(res) != 0) {
         begin <- end(tail(res, n=-1))
      }
      if (i %in% subjectHits(hits_int)) {
         hit <- hits_int[subjectHits(hits_int)==i]
         if (queryHits(hit) == 1) {
            res <- c(res, GRanges(seqnames=chr, strand=strand, 
                                  ranges=IRanges(start=start(subjectHits(hit)),
                                                 end=end(queryHits(hit))),
                                  type='a')
         } else if (queryHits(hit) == last_int) {
            res <- c(res, GRanges(seqnames=chr, strand=strand, 
                                  ranges=IRanges(start=start(queryHits(hit)),
                                                 end=end(subjectHits(hit))),
                                  type='a')
         } else {
            res <- c(res, GRanges(seqnames=chr, strand=strand, 
                                  ranges=IRanges(start=begin),
                                                 end=start(queryHits(hit))),
                                  type='I')
            res <- c(res, GRanges(seqnames=chr, strand=strand, 
                                  ranges=IRanges(start=end(queryHits(hit)),
                                                 end=end(subjectHits(hit))),
                                  type='A')
         }
      } else if (i %in% subjectHits(hits_ex)) {
         hit <- hits_ex[subjectHits(hits_ex)==i]
         type <- 'e'
      }
   }
   
}