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
   strand <- unique(as.character(strand(geneGr)))
   chr <- unique(as.character(seqnames(geneGr)))
   gene_id <- unique(as.character(mcols(geneGr)$gene))
   if (length(chr) != 1) {
      stop("A gene is on different chrs")
   }
   if (length(gene_id) != 1) {
      stop("A gene has different names")
   }
   if (strand == '+') {
      getApaGenesFractionsPlusStrand(geneGr, apaGr, chr, strand, gene_id)
   } else if (strand == '-') {
      stop("Implementation for minus strand is still missing")
      # TODO minus varsion
   } else {
      stop("A gene has no strand info or both + and -")
   } 
}

# XXX we need to have the seqlengths to have a correct last intron. TODO
# maybe we could avoid doing this...we need to have faith in the gtf tough
getApaGenesFractionsPlusStrand <- function(geneGr, apaGr, chr, strand, gene_id)
{   
   geneGr <- sort(geneGr)
   apaGr <- sort(apaGr)
   introns <- gaps(geneGr)
   # Remove first intron: it is not useful in any way.
   introns <- tail(introns, n=length(introns)-1)
   lastApa <- tail(apaGr, n=1)
   names <- unlist(strsplit(mcols(lastApa)$apa, '_', fixed=TRUE))
   lastApa_name <- names[1]
   apaFragmentsPrePost <- vector("list", length(apaGr)-1)
   apaFrI <- 1
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how. 
   # New idea: we combine exons and introns adding them mcols to their kind
   # and with info about being overlappant or not.
   # Then we run over this GRanges that should have all the needed info
   # to generate a GRanges res with intervals and maybe also a list of 
   # pre/post indexes for every APA?
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'   
   if (length(introns) != 0) {
      mcols(introns)$type <- 'i'  
   }   
   whole <- sort(c(geneGr, introns))
   hits <- findOverlaps(whole, apaGr)
   # We have downstream apas so we can't do this.
   #if (!all(countSubjectHits(hits) == 1)) {
   #   stop("Error: a given apa does not overlap its gene")
   #}
   
   foundov <- length(unique(whole[queryHits(hits)]))
   mcols(whole)$overlap <- rep(FALSE, length(whole))   
   mcols(whole[unique(queryHits(hits))])$overlap <- rep(TRUE, foundov)
   # Growing objects sin.
   begins <- c()
   ends <- c()
   begin <- -1
   lastExonB <- -1
   # Write logic and invariant here: begin is always the beginning of the last
   # needed fragment (i.e. start of exon with an overlap or with an overlap in
   # the next intron or just seen apa). We add an end when 
   # we get an exon/intron with overlapping apas, and a begin for exons
   # with overlap or preceding an overlapping intron.
   for (i in 1:length(whole)) {
      if (whole[i]$type=='e') {
         if (whole[i]$overlap || (i+1<=length(whole) && whole[i+1]$overlap)) {
            if (begin != -1) {
               begins <- c(begins, begin)
               ends <- c(ends, end(whole[i-1]))
            }
            begin <- start(whole[i])
            lastExonB <- length(begins)+1
         }
      }
      if (whole[i]$overlap) {
         overlapping_apas <- apaGr[subjectHits(hits[queryHits(hits)==i])]
         # Here we have a shot of writing the pre portions of counts
         # for each APA.
         for (j in 1:length(overlapping_apas)) {
            begins <- c(begins, begin)
            ends <- c(ends, start(overlapping_apas[j]))
            # We add 1 to avoid overlapping PRE/POST fragments 
            # - APA are considered before the cut.
            # In the previous gtf we skipped the "cut" base altogether.
            # XXX need to decide.
            begin <- end(overlapping_apas[j])+1
            firstApa_name <- unlist(
               strsplit(mcols(overlapping_apas[j])$apa, '_', fixed=TRUE))[1]
            apa_name <- paste(firstApa_name, lastApa_name, sep="-")
            apaFragmentsPrePost[[apaFrI]] <- new("ApaFragmentPrePost",
                                                   name=apa_name,
                                                   PREstart=lastExonB, 
                                                   PREend=length(begins))
            apaFrI <- apaFrI+1
         }
      }
   }
   # We could be ok like this (last APA inside an exon, we stop there)
   # or we could have a last APA after the last exon that we want to consider.
   #if (an APA outside hits) {
   if (!all(countSubjectHits(hits) == 1)) {
      out_apas <- apaGr[countSubjectHits(hits) == 0]
      # If we did not see an overlap (should be rare) we have a single fragment 
      # starting at the beginning of the last exon.
      if (begin == -1) {
         begin <- start(tail(geneGr, n=1)) 
         lastExonB <- length(begins)+1 # That's always 1 I know.
      } 
      for (j in 1:length(out_apas)) {
         begins <- c(begins, begin)
         ends <- c(ends, start(out_apas[j]))
         # We add 1 to avoid overlapping PRE/POST fragments 
         # - APA are considered before the cut.
         # In the previous gtf we skipped the "cut" base altogether.
         # XXX need to decide.
         begin <- end(out_apas[j])+1
         firstApa_name <- unlist(
            strsplit(mcols(out_apas[j])$apa, '_', fixed=TRUE))[1]
         apa_name <- paste(firstApa_name, lastApa_name, sep="-")
         apaFragmentsPrePost[[apaFrI]] <- new("ApaFragmentPrePost",
                                              name=apa_name,
                                              PREstart=lastExonB, 
                                              PREend=length(begins))
         apaFrI <- apaFrI+1
      }
   }
   fragments <- GRanges(seqnames=chr, strand=strand, 
                        ranges=IRanges(start=begins, end=ends))
   # The last apaFragment is build from a single apa only and we don't want it
   #fragments <- head(fragments, n=length(fragments)-1)
   # But the fragment is ok!
   apaFragmentsPrePost <- head(apaFragmentsPrePost,
                              n=apaFrI-2)
   return(list(fragments, apaFragmentsPrePost))
}
# Testing will be hideous.

obtainPrePost <- function(prepost, fragments)
{
   sn <- unique(seqnames(fragments))
   gene_id=c(paste0(prepost@name, "_PRE"),
             paste0(prepost@name, "_POST"))
   # We do not really need the right coords, but still: I was fool.
   # We need them to correct in the right way for lengths!
   res <- GRanges(seqnames=rep(sn,2),
         ranges=IRanges(start=c(start(fragments[prepost@PREstart]),
                                start(fragments[prepost@PREend+1])),
                        end=c(end(fragments[prepost@PREend]),
                              end(tail(fragments, n=1)))))
   mcols(res) <- DataFrame(gene_id)
   return(res)
   # XXX TODO strand minus.
}

sumFragmentCounts<- function(prepost, counts, kind)
{
   if (kind == "pre") 
   {
      return(sum(counts[c(prepost@PREstart, prepost@PREend)]))
   } else if (kind == "post") {
      return(sum(counts[c(prepost@PREend+1, length(counts))]))
   } else {
      stop ("Inner error in helper functions to sum counts")
   }
}

createRoarSingleBam <- function(name, mulRds, treatmentSE, controlSE)
{
   prePostDef <- mulRds@prePostDef[[name]]
   fragments <- mulRds@fragments[[name]]
   #prePostCoords = "GRanges",
   #postCoords = "GRanges",
   #countsTreatment = "RangedSummarizedExperiment",
   #countsControl = "RangedSummarizedExperiment",
   prePostCoords <- do.call(c, sapply(prePostDef, obtainPrePost, fragments))
   # TODO XXX add entrez_id here! We have a roar object foreach gene
   # so maybe we do not need them.
   # But we need to extract the right values from treatmentSE/controlSE, should
   # clearly be more efficient than looking for the coords in fragments.
   #  assays(treatmentSE,1)$counts[rownames(assays(treatmentSE,1)$counts)=="10771"]
   postElems <- grep("_POST$", mcols(prePostCoords)$gene_id)
   preElems <- grep("_PRE$", mcols(prePostCoords)$gene_id)
   preCoords <- prePostCoords[preElems,]
   rds <- new("RoarDataset", prePostCoords=prePostCoords)
   rds@postCoords <- rds@prePostCoords[postElems,]
   se <- SummarizedExperiment(assays = rep(list(matrix(nrow=length(rds@prePostCoords)/2, ncol=4)),2),
                              rowRanges=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   treatment <- assays(treatmentSE,1)$counts[rownames(assays(treatmentSE,1)$counts)==name]
   control <- assays(controlSE,1)$counts[rownames(assays(controlSE,1)$counts)==name]
   assay(se,1)[,"treatment_pre"] <- sapply(prePostDef, sumFragmentCounts, treatment, "pre")
   assay(se,1)[,"treatment_post"] <- sapply(prePostDef, sumFragmentCounts, treatment, "post")
   assay(se,1)[,"control_pre"] <- sapply(prePostDef, sumFragmentCounts, control, "pre")
   assay(se,1)[,"control_post"] <- sapply(prePostDef, sumFragmentCounts, control, "post")
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@step <- 1
   rds@treatmentBams <- mulRds@treatmentBams
   rds@controlBams <- mulRds@controlBams
   rds@cores <- 1
   return(rds)
}