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
      getApaGenesFractionsMinusStrand(geneGr, apaGr, chr, strand, gene_id)
   } else {
      stop("A gene has no strand info or both + and -")
   } 
}

getApaGenesFractionsPlusStrand <- function(geneGr, apaGr, chr, strand, gene_id)
{   
   geneGr <- sort(geneGr)
   apaGr <- sort(apaGr)
   seqlengths(seqinfo(geneGr))<-rep(NA, length(seqinfo(geneGr)))
   introns <- gaps(geneGr)
   # Remove first intron: it is not useful in any way.
   # TODO: if there are chr lenghts we will get also a last intron? Yes. Bad fix: remove them.
   introns <- tail(introns, n=length(introns)-1) # TODO add check to avoid removing introns if first exon begin == 0/1
   fakeApaEndGene <- GRanges(
      seqnames = chr,
      strand = strand,
      ranges = IRanges(
         start = end(tail(geneGr, n=1)),
         width = 1,
         names = ""
      )
   )
   mcols(fakeApaEndGene) <- mcols(apaGr)[1,]
   mcols(fakeApaEndGene)[1,"apa"] <- ""
   apaGr <- c(apaGr, fakeApaEndGene)
   apaFragmentsPrePost <- vector("list", length(apaGr)-1)
   apaFrI <- 1
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how. 
   # New idea: we combine exons and introns adding them mcols to their kind
   # and with info about being overlappant or not.
   # Then we run over this GRanges that should have all the needed info
   # to generate a GRanges res with intervals and maybe also a list of 
   # pre/post indexes for every APA.
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'
   if (length(introns) != 0) {
      mcols(introns)$type <- 'i'  
   }   
   whole <- sort(c(geneGr, introns))
   hits <- findOverlaps(whole, apaGr)
   if (!all(countSubjectHits(hits) == 1)) {
      stop("Error: a given apa does not overlap its gene or an apa is not of width 1")
   }
   
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
            begin <- end(overlapping_apas[j])+1
            firstApa_name <- unlist(
               strsplit(mcols(overlapping_apas[j])$apa, '_', fixed=TRUE))[1]
            apa_name <- firstApa_name
            apaFragmentsPrePost[[apaFrI]] <- new("ApaFragmentPrePost",
                                                   name=apa_name,
                                                   PREstart=lastExonB, 
                                                   PREend=length(begins))
            apaFrI <- apaFrI+1
         }
      }
   }
   fragments <- GRanges(seqnames=chr, strand=strand, 
                        ranges=IRanges(start=begins, end=ends))
   ovHits <- findOverlaps(fragments, geneGr)
   ovLen <- width(overlapsRanges(ranges(fragments), ranges(geneGr), hits=ovHits))
   mcols(fragments)$length <- rep(0, length(fragments))
   # A single fragment can overlap multiple exons and we need to sum here their lengths:
   # i.e. ovlen is longer than length(fragments), with queryHits we get the right reps but
   # not sums. The only way I am able to imagine a solution is with a for.
   indexes_fr <- queryHits(ovHits)
   for (k in 1:length(ovLen)) {
      mcols(fragments[indexes_fr[k]])$length  <- mcols(fragments[indexes_fr[k]])$length + ovLen[k]
   }
   # The last apaFragment is build from a single apa (end of gene) only and we don't want it
   #fragments <- head(fragments, n=length(fragments)-1)
   # But the fragment is ok!
   apaFragmentsPrePost <- head(apaFragmentsPrePost,
                              n=apaFrI-2)
   return(list(fragments, apaFragmentsPrePost))
}
# Testing will be hideous.

# Strand -. We are being not DRY but i prefer to start in this way and refactor
# after the first runs/tests. During active development I will have to remember
# to change code (if needed) in two places...but right now I'm not sure
# about how many differences I need to insert in the code.
getApaGenesFractionsMinusStrand <- function(geneGr, apaGr, chr, strand, gene_id)
{   
   geneGr <- sort(geneGr, decreasing=TRUE)
   apaGr <- sort(apaGr, decreasing=TRUE)
   introns <- sort(gaps(geneGr), decreasing=TRUE)
   seqlengths(seqinfo(geneGr))<-rep(NA, length(seqinfo(geneGr)))
   # Remove last intron: it is not useful in any way.
   introns <- head(introns, n=length(introns)-1)
   fakeApaEndGene <- GRanges(
      seqnames = Rle(c(chr)),
      strand = strand,
      ranges = IRanges(
         start= start(tail(geneGr, n=1)),
         width= 1,
         names = ""
      )
   )
   mcols(fakeApaEndGene) <- mcols(apaGr)[1,]
   mcols(fakeApaEndGene)[1,"apa"] <- ""
   apaGr <- c(apaGr, fakeApaEndGene)
   apaFragmentsPrePost <- vector("list", length(apaGr)-1)
   apaFrI <- 1
   # We start by doing it iteratively and non R/Bioc stylishly because I've got
   # no idea how. 
   # New idea: we combine exons and introns adding them mcols to their kind
   # and with info about being overlappant or not.
   # Then we run over this GRanges that should have all the needed info
   # to generate a GRanges res with intervals and maybe also a list of 
   # pre/post indexes for every APA.
   mcols(geneGr) <- NULL
   mcols(geneGr)$type <- 'e'   
   if (length(introns) != 0) {
      mcols(introns)$type <- 'i'  
   }   
   whole <- sort(c(geneGr, introns), decreasing=TRUE)
   hits <- findOverlaps(whole, apaGr)
   if (!all(countSubjectHits(hits) == 1)) {
      stop("Error: a given apa does not overlap its gene or an apa is not of width 1")
   }
   
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
               ends <- c(ends, start(whole[i-1]))
            }
            begin <- end(whole[i])
            lastExonB <- length(begins)+1
         }
      }
      if (whole[i]$overlap) {
         overlapping_apas <- apaGr[subjectHits(hits[queryHits(hits)==i])]
         overlapping_apas <- sort(overlapping_apas, decreasing=TRUE)
         # Here we have a shot of writing the pre portions of counts
         # for each APA.
         for (j in 1:length(overlapping_apas)) {
            begins <- c(begins, begin)
            ends <- c(ends, start(overlapping_apas[j]))
            # We add 1 to avoid overlapping PRE/POST fragments 
            # - APA are considered before the cut.
            # In the previous gtf we skipped the "cut" base altogether.
            begin <- end(overlapping_apas[j])-1
            firstApa_name <- unlist(
               strsplit(mcols(overlapping_apas[j])$apa, '_', fixed=TRUE))[1]
            apa_name <- firstApa_name
            apaFragmentsPrePost[[apaFrI]] <- new("ApaFragmentPrePost",
                                                 name=apa_name,
                                                 PREstart=lastExonB, 
                                                 PREend=length(begins))
            apaFrI <- apaFrI+1
         }
      }
   }
   # here begins and ends are to be switched.
   fragments <- GRanges(seqnames=chr, strand=strand, 
                        ranges=IRanges(start=ends, end=begins))
   ovHits <- findOverlaps(fragments, geneGr)
   ovLen <- width(overlapsRanges(ranges(fragments), ranges(geneGr), hits=ovHits))
   
   mcols(fragments)$length <- rep(0, length(fragments))
   
   # A single fragment can overlap multiple exons and we need to sum here their lengths:
   # i.e. ovlen is longer than length(fragments), with queryHits we get the right reps but
   # not sums. The only way I am able to imagine a solution is with a for.
   indexes_fr <- queryHits(ovHits)
   for (k in 1:length(ovLen)) {
      mcols(fragments[indexes_fr[k]])$length  <- mcols(fragments[indexes_fr[k]])$length + ovLen[k]
   }
   #mcols(fragments[queryHits(ovHits)])$length  <- ovLen # But for length is here XXX FIXME
   
   # The last apaFragment is build from a single apa only and we don't want it
   #fragments <- head(fragments, n=length(fragments)-1)
   # But the fragment is ok!
   apaFragmentsPrePost <- head(apaFragmentsPrePost,
                               n=apaFrI-2)
   return(list(fragments, apaFragmentsPrePost))
}


obtainPrePost <- function(prepost, fragments)
{
   sn <- unique(seqnames(fragments))
   strand <- unique(strand(fragments)) 
   # We have already checked about strand uniqueness.
   gene_id=c(paste0(prepost@name, "_PRE"),
             paste0(prepost@name, "_POST"))
   # We do not really need the right coords, but still: I was fool.
   # We need them to correct in the right way for lengths!
   if (strand == "+") {
      # This logic does not work for minus strand.
      res <- GRanges(seqnames=rep(sn,2), # 22% of the time here?
                     strand=rep(strand,2),
            ranges=IRanges(start=c(start(fragments[prepost@PREstart]),
                                   start(fragments[prepost@PREend+1])),
                           end=c(end(fragments[prepost@PREend]),
                                 end(tail(fragments, n=1)))))
      length <- c(sum(mcols(fragments[seq(prepost@PREstart, prepost@PREend)])$length),
                             sum(mcols(fragments[seq(prepost@PREend+1, length(fragments))])$length))
   } else { # if (strand == "-") { # and also only + and -
      res <- GRanges(seqnames=rep(sn,2), # 22% of the time here?
                     strand=rep(strand,2),
                     ranges=IRanges(start=c(start(fragments[prepost@PREend]),
                                            start(tail(fragments, n=1))),
                                    end=c(end(fragments[prepost@PREstart]),
                                          end(fragments[prepost@PREend+1]))))
      length <- c(sum(mcols(fragments[seq(prepost@PREstart, prepost@PREend)])$length),
                  sum(mcols(fragments[seq(prepost@PREend+1, length(fragments))])$length))
   }
   mcols(res) <- DataFrame(gene_id, length)
   return(res)
}

sumFragmentCounts<- function(prepost, counts, kind)
{
   if (kind == "pre") {
      return(sum(counts[c(seq(prepost@PREstart,prepost@PREend, by=1))]))
   } else if (kind == "post") {
      return(sum(counts[c(seq(prepost@PREend+1,length(counts), by=1))]))
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
   # We need to extract the right values from treatmentSE/controlSE, should
   # clearly be more efficient than looking for the coords in fragments.
   #  assays(treatmentSE,1)$counts[rownames(assays(treatmentSE,1)$counts)=="10771"]
   postElems <- grep("_POST$", mcols(prePostCoords)$gene_id)
   preElems <- grep("_PRE$", mcols(prePostCoords)$gene_id)
   preCoords <- prePostCoords[preElems,]
   rds <- new("RoarDataset", 
              SummarizedExperiment(rowRanges=preCoords,
                                   colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))),
              prePostCoords=prePostCoords)
   rds@postCoords <- rds@prePostCoords[postElems,]
   se <- SummarizedExperiment(assays = rep(list(matrix(nrow=length(rds@prePostCoords)/2, ncol=4)),2),
                              rowRanges=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   # Other than assay rownames we could have used:
   #> rownames(assays(pada,1)$counts)
   #[1] "100996928" "100996928" "101927572" "101927572" "104909134" "104909134"
   #[7] "9"         "9"         "10"        "10"        "12"        "12"       
   #> names(rowRanges(pada))
   #[1] "100996928" "100996928" "101927572" "101927572" "104909134" "104909134"
   #[7] "9"         "9"         "10"        "10"        "12"        "12"      
   #treatment <- assays(treatmentSE,1)$counts[rownames(assays(treatmentSE,1)$counts)==name]
   # > assays(pada,1)$counts[names(rowRanges(pada))=="12"]
   # [1] 0 0
   # But which is faster?
   control <- assays(controlSE,1,withDimnames=TRUE)$counts[rownames(assays(controlSE,1,withDimnames=TRUE)$counts)==name]
   treatment <- assays(treatmentSE,1,withDimnames=TRUE)$counts[rownames(assays(treatmentSE,1,withDimnames=TRUE)$counts)==name]
   assay(se,1, withDimnames=TRUE)[,"treatment_pre"] <- sapply(prePostDef, sumFragmentCounts, treatment, "pre")
   assay(se,1, withDimnames=TRUE)[,"treatment_post"] <- sapply(prePostDef, sumFragmentCounts, treatment, "post")
   assay(se,1, withDimnames=TRUE)[,"control_pre"] <- sapply(prePostDef, sumFragmentCounts, control, "pre")
   assay(se,1, withDimnames=TRUE)[,"control_post"] <- sapply(prePostDef, sumFragmentCounts, control, "post")
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds, withDimnames=TRUE) <- assays(se, withDimnames=TRUE)
   names(assays(rds,withDimnames=TRUE)) <- "counts"
   rds@step <- 1
   # length should be preserved: XXX multipleSamples!
   rds@treatmentBams <- list(NA)
   rds@controlBams <- list(NA)
   rds@cores <- 1
   return(rds)
}

createRoarMultipleBam <- function(name, mulRds, treatmentSE, controlSE)
{
   prePostDef <- mulRds@prePostDef[[name]]
   fragments <- mulRds@fragments[[name]]
   #prePostCoords = "GRanges",
   #postCoords = "GRanges",
   #countsTreatment = "RangedSummarizedExperiment",
   #countsControl = "RangedSummarizedExperiment",
   prePostCoords <- do.call(c, sapply(prePostDef, obtainPrePost, fragments))
   # We need to extract the right values from treatmentSE/controlSE, should
   # clearly be more efficient than looking for the coords in fragments.
   #  assays(treatmentSE,1)$counts[rownames(assays(treatmentSE,1)$counts)=="10771"]
   postElems <- grep("_POST$", mcols(prePostCoords)$gene_id)
   preElems <- grep("_PRE$", mcols(prePostCoords)$gene_id)
   preCoords <- prePostCoords[preElems,]
   rds <- new("RoarDataset", 
              SummarizedExperiment(rowRanges=preCoords,
                                   colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))),
              prePostCoords=prePostCoords)
   rds@postCoords <- rds@prePostCoords[postElems,]
   se <- SummarizedExperiment(assays = rep(list(matrix(nrow=length(rds@prePostCoords)/2, ncol=4)),2),
                              rowRanges=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rds@countsTreatment <- SummarizedExperiment(assays = rep(list(matrix(nrow=length(rds@prePostCoords)/2, ncol=2)), length(mulRds@treatmentBams)),
                              rowRanges=preCoords, 
                              colData=DataFrame(row.names=c("pre","post"))
   )
   rds@countsControl <- SummarizedExperiment(assays = rep(list(matrix(nrow=length(rds@prePostCoords)/2, ncol=2)), length(mulRds@controlBams)),
                                           rowRanges=preCoords, 
                                           colData=DataFrame(row.names=c("pre", "post"))
   )
   
   for (i in 1:length(controlSE)) {
      control <- assays(controlSE[[i]],1,withDimnames=TRUE)$counts[rownames(assays(controlSE[[i]],1,withDimnames=TRUE)$counts)==name]
      assay(rds@countsControl,i,withDimnames=TRUE)[,"pre"] <- sapply(prePostDef, sumFragmentCounts, control, "pre")
      assay(rds@countsControl,i,withDimnames=TRUE)[,"post"] <- sapply(prePostDef, sumFragmentCounts, control, "post")
   }
   for (i in 1:length(treatmentSE)) {
      treatment <- assays(treatmentSE[[i]],1,withDimnames=TRUE)$counts[rownames(assays(treatmentSE[[i]],1,withDimnames=TRUE)$counts)==name]
      assay(rds@countsTreatment,i,withDimnames=TRUE)[,"pre"] <- sapply(prePostDef, sumFragmentCounts, treatment, "pre")
      assay(rds@countsTreatment,i,withDimnames=TRUE)[,"post"] <- sapply(prePostDef, sumFragmentCounts, treatment, "post")
   }
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds,withDimnames=TRUE) <- assays(se,withDimnames=TRUE)
   names(assays(rds,withDimnames=TRUE)) <- "counts"
   rds@step <- 1
   rds@treatmentBams <- vector(mode = "list", length = length(mulRds@treatmentBams))
   rds@controlBams <- vector(mode = "list", length = length(mulRds@controlBams))   
   rds@cores <- 1
   return(rds)
}
