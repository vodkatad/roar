# To manually perform single tests: source("inst/unitTests/test_roar.R")
# To perform all tests: 
# library("GenomicAlignments")
# BiocGenerics:::testPackage("roar")
# countPrePost ------------------------------------------------------------

test_countPrePost_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", "D_POST", "E_PRE", "E_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2", "chr2", "chr2",
                       "chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000),
         width=c(500, 900, 500, 300, 600, 300, 500, 900)),
      DataFrame(gene_id)
   )
   rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(c(rd1,rd2)), list(c(rd3,rd4)), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds,1)[1,1], 1)
   checkEqualsNumeric(assay(rds,1)[1,2], 1)
   checkEqualsNumeric(assay(rds,1)[2,3], 2)
}

test_countPrePost_mulSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", "D_POST", "E_PRE", "E_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2", "chr2", "chr2",
                       "chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000),
         width=c(500, 900, 500, 300, 600, 300, 500, 900)),
      DataFrame(gene_id)
   )
   rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(rd1,rd2), list(rd3,rd4), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds@countsTreatment,1)[1,1], 1)
   checkEqualsNumeric(assay(rds@countsTreatment,1)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsControl,2)[2,1], 1)
}

test_countPrePost_preferPOST <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   # Had to change from Rle() to factor for seqnames otherwise there was a warning about the list for the control
   # alignments about having only chr1 for some reads and chr2 for the others.
   rd1 <- GAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
   rd2 <- GAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
   rd3 <- GAlignments("a", seqnames = factor("chr2", levels=c("chr1","chr2")), pos = as.integer(2800), cigar = "300M", strand = strand("+"))
   rd4 <- GAlignments("a", seqnames = factor("chr2", levels=c("chr1","chr2")), pos = as.integer(3500), cigar = "300M", strand = strand("+"))
   rd5 <- GAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(c(rd1,rd2)), list(c(rd3,rd4,rd5)), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds,1)[1,1], 0)
   checkEqualsNumeric(assay(rds,1)[1,2], 2)
   checkEqualsNumeric(assay(rds,1)[1,3], 1)
   checkEqualsNumeric(assay(rds,1)[1,4], 0)
   checkEqualsNumeric(assay(rds,1)[2,3], 1)
   checkEqualsNumeric(assay(rds,1)[2,4], 1)
   checkEqualsNumeric(assay(rds,1)[2,1], 0)
   checkEqualsNumeric(assay(rds,1)[2,2], 0)
}

test_countPrePost_stranded <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("-"))
   rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("-"))
   rd5 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   rd6 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(rd1,rd2), list(rd3,rd4,c(rd5,rd6)), features)
   rds <- countPrePost(rds, TRUE)
   checkEqualsNumeric(assay(rds@countsTreatment,1)[1,1], 1)
   checkEqualsNumeric(assay(rds@countsTreatment,1)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsTreatment,2)[1,1], 0)
   checkEqualsNumeric(assay(rds@countsTreatment,2)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsControl,1)[2,1], 1)
   checkEqualsNumeric(assay(rds@countsControl,2)[2,1], 0)
   checkEqualsNumeric(assay(rds@countsControl,3)[2,1], 2)
}

# These UnitTests are bad because I should "prebuild" a correct step 1 complete
# rds object to test computeRoars and avoid testing multiple features at the same time,
# but it is too complex and overall the tests should work even in this way.

# computeRoars ------------------------------------------------------------
test_computeRoars_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_POST", "B_PRE", "C_PRE", "C_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(c("+","+","-","-","+","+")),
      ranges = IRanges(
         start=c(1, 10, 20, 40, 42, 52),
         width=c(10, 5, 20, 5, 10, 10)),
      DataFrame(gene_id)
   )
   a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), cigar = "5M", strand = strand("+"))
   a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), cigar = "3M", strand = strand("+"))
   a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), cigar = "3M", strand = strand("+"))
   b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "1M", strand = strand("-"))
   b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), cigar = "3M", strand = strand("-"))
   # The next one is an overlapping read only if strandness is not considered. 
   # Otherwise it will be counted for pre_C. I will add a check on this in this test, even if it's not correct.
   overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "5M", strand = strand("+"))
   c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), cigar = "50M", strand = strand("+"))
   
   
   treatmentAlign <- list(c(rep(a_pre, 2),rep(a_post, 3), a_pre_post, rep(b_pre,5), b_post, overlapbc))
   controlAlign <- list(c(a_post, rep(a_pre, 4), a_pre_post, rep(b_post,5), b_pre, c_post))
   
   rds <- RoarDataset(treatmentAlign, controlAlign, features)
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
   #assay(rds,2) <- as.matrix(data.frame(treatment_pre=mMtreatment, treatment_post=mMcontrol, control_pre=roar, control_post=pVal))
   checkEqualsNumeric(assay(rds,2)[1,1], -0.66538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,2], 1.21538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,3], -0.54746835443038, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[1,4]))
   checkEqualsNumeric(assay(rds,2)[2,1], 20.6923076923077, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,2], 0.04307692307692, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,3], 480.357142857144, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[3,1]))
   checkEqualsNumeric(assay(rds,2)[3,2], -1)
   checkTrue(is.na(assay(rds,2)[3,3]))      
   
   # The out of place test about stranded alignment:
   rds2 <- RoarDataset(treatmentAlign, controlAlign, features)
   rds2 <- countPrePost(rds2, TRUE)
   checkEqualsNumeric(assay(rds2,1)[1,1],2)
   checkEqualsNumeric(assay(rds2,1)[1,2],4)
   checkEqualsNumeric(assay(rds2,1)[1,3],4)
   checkEqualsNumeric(assay(rds2,1)[1,4],2)
   checkEqualsNumeric(assay(rds2,1)[2,1],5)
   checkEqualsNumeric(assay(rds2,1)[2,2],1)
   checkEqualsNumeric(assay(rds2,1)[2,3],1)
   checkEqualsNumeric(assay(rds2,1)[2,4],5)
   checkEqualsNumeric(assay(rds2,1)[3,1],1)
   checkEqualsNumeric(assay(rds2,1)[3,2],0)
   checkEqualsNumeric(assay(rds2,1)[3,3],0)
   checkEqualsNumeric(assay(rds2,1)[3,4],1)
}

test_computeRoars_singlevsMulSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_POST", "B_PRE", "C_PRE", "C_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(c("+","+","-","-","+","+")),
      ranges = IRanges(
         start=c(1, 10, 20, 40, 42, 52),
         width=c(10, 5, 20, 5, 10, 10)),
      DataFrame(gene_id)
   )
   a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), cigar = "5M", strand = strand("+"))
   a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), cigar = "3M", strand = strand("+"))
   a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), cigar = "3M", strand = strand("+"))
   b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "1M", strand = strand("-"))
   b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), cigar = "3M", strand = strand("-"))
   overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "5M", strand = strand("+"))
   c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), cigar = "50M", strand = strand("+"))
   
   treatmentAlign <- list(c(rep(a_pre, 2),rep(a_post, 3), a_pre_post, rep(b_pre,5), b_post, overlapbc))
   controlAlign <- list(a_post, rep(a_pre, 4), a_pre_post, rep(b_post,5), b_pre, c_post)
   
   rds <- RoarDataset(treatmentAlign, controlAlign, features)
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
   #assay(rds,2) <- as.matrix(data.frame(treatment_pre=mMtreatment, treatment_post=mMcontrol, control_pre=roar, control_post=pVal))
   checkEqualsNumeric(assay(rds,2)[1,1], -0.66538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,2], 1.21538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,3], -0.54746835443038, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[1,4]))
   checkEqualsNumeric(assay(rds,2)[2,1], 20.6923076923077, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,2], 0.04307692307692, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,3], 480.357142857144, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[3,1]))
   checkEqualsNumeric(assay(rds,2)[3,2], -1)
   checkTrue(is.na(assay(rds,2)[3,3]))      
}

test_computeRoars_multipleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1","chr1")),
      strand = strand(c("+","+")),
      ranges = IRanges(
         start=c(1, 10),
         width=c(10, 5)),
      DataFrame(gene_id)
   )
   a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), cigar = "5M", strand = strand("+"))
   a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), cigar = "3M", strand = strand("+"))
      
   treatmentAlign <- list(rep(a_pre, 2),rep(a_post, 3))
   controlAlign <- list(a_post, rep(a_pre, 4), a_post)
   
   rds <- RoarDataset(treatmentAlign, controlAlign, features)
   rds <- countPrePost(rds, FALSE)
   rds <- computeRoars(rds)
   #assay(rds,2) <- as.matrix(data.frame(treatment_pre=mMtreatment, treatment_post=mMcontrol, control_pre=roar, control_post=pVal))
   checkEqualsNumeric(assay(rds,2)[1,1], -0.48, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,2], 0.666666666, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,3], -0.72, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[1,4]))
}

# computePvals ------------------------------------------------------------
# Ok, to test computePVals I know how to fill the rds object with counts so I will skip the previous
# steps and perform slightly more correct tests.

test_computePvals_single <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
       prePostCoords=features, step = 2, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   # We need to set these lengths to choose the treatment branch of the if in computePvals (about number of samples).
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 1
   assay(rds, 1)[1,1] <- 10
   assay(rds, 1)[1,2] <- 100
   assay(rds, 1)[1,3] <- 20
   assay(rds, 1)[1,4] <- 20
   assay(rds, 1)[2,1] <- 10
   assay(rds, 1)[2,2] <- 20
   assay(rds, 1)[2,3] <- 20
   assay(rds, 1)[2,4] <- 20
   # We need to setup the second assay that computePVals will fill.
   assay(rds,2) <- matrix(nrow=2, ncol=4)
   rds <- computePvals(rds)
   checkEqualsNumeric(assay(rds,2)[1,4], 0.0000002212406, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,4], 0.22337243, tolerance=1e-5)
}

test_computePvals_singlevsMul <- function() {
   gene_id <- c("A_PRE", "A_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000),
         width=c(1000, 900)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 2, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=1, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   # We need to set these lengths to choose the treatment branch of the if in computePvals (about number of samples).
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 3
   rds@countsTreatment <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                                        rowData=preCoords, 
                                                        colData=DataFrame(row.names=c("pre","post"))
   )
   rds@countsControl <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                       rowData=preCoords, 
                                       colData=DataFrame(row.names=c("pre","post"))
   )
   assay(rds@countsControl,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsControl,3) <- matrix(nrow=1,ncol=2)
   assay(rds@countsTreatment,1)[1,] <- c(10,10)
   assay(rds@countsControl,1)[1,] <- c(10,10)
   assay(rds@countsControl,2)[1,] <- c(10,5)
   assay(rds@countsControl,3)[1,] <- c(10,20)
   # We need to setup the second assay that computePVals will fill.
   assay(rds,2) <- matrix(nrow=1, ncol=4)
   rds <- computePvals(rds)
   checkEqualsNumeric(assay(rds@pVals,1)[1,1], 1, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,2], 0.491599, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,3], 0.257549, tolerance=1e-5)
}

test_computePvals_multipleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000),
         width=c(1000, 900)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 2, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=1, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   # We need to set these lengths to choose the right branch of the if in computePvals (about number of samples).
   length(rds@treatmentBams)  <- 3
   length(rds@controlBams)  <- 2
   rds@countsTreatment <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                           rowData=preCoords, 
                                           colData=DataFrame(row.names=c("pre","post"))
   )
   rds@countsControl <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                          rowData=preCoords, 
                                          colData=DataFrame(row.names=c("pre","post"))
   )
   assay(rds@countsTreatment,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsTreatment,3) <- matrix(nrow=1,ncol=2)
   assay(rds@countsControl,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsTreatment,1)[1,] <- c(10,10)
   assay(rds@countsTreatment,2)[1,] <- c(1,10)
   assay(rds@countsTreatment,3)[1,] <- c(20,10)
   assay(rds@countsControl,1)[1,] <- c(10,5)
   assay(rds@countsControl,2)[1,] <- c(10,20)
   # 6 pvalues: 1_1 -> 0.491599
   # 1_2 -> 0.257549
   # 2_1 -> 0.00522109
   # 2_2 -> 0.23264154
   # 3_1 -> 1
   # 3_2 -> 0.01938319
   # We need to setup the second assay that computePVals will fill.
   assay(rds,2) <- matrix(nrow=1, ncol=4)
   rds <- computePvals(rds)
   checkEqualsNumeric(assay(rds@pVals,1)[1,1], 0.491599, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,2], 0.257549, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,3], 0.00522109, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,4], 0.23264154, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,5], 1, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,6], 0.01938319, tolerance=1e-5)
}

test_computePairedPvals <- function() {
   gene_id <- c("A_PRE", "A_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000),
         width=c(1000, 900)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 2, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=1, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   
   rowRanges(rds) <- rowRanges(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   # We need to set these lengths to choose the treatment branch of the if in computePvals (about number of samples).
   length(rds@treatmentBams)  <- 3
   length(rds@controlBams)  <- 2
   rds@countsTreatment <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                               rowData=preCoords, 
                                               colData=DataFrame(row.names=c("pre","post"))
   )
   rds@countsControl <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                             rowData=preCoords, 
                                             colData=DataFrame(row.names=c("pre","post"))
   )
   assay(rds@countsTreatment,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsTreatment,3) <- matrix(nrow=1,ncol=2)
   assay(rds@countsControl,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsTreatment,1)[1,] <- c(10,10)
   assay(rds@countsTreatment,2)[1,] <- c(1,10)
   assay(rds@countsTreatment,3)[1,] <- c(20,10)
   assay(rds@countsControl,1)[1,] <- c(10,5)
   assay(rds@countsControl,2)[1,] <- c(10,20)
   # 2_1 -> 0.00522109
   # 3_2 -> 0.01938319
   # We need to setup the second assay that computePairedVals will fill.
   assay(rds,2) <- matrix(nrow=1, ncol=4)
   rds <- computePairedPvals(rds, c(2,3), c(1,2))
   checkEqualsNumeric(assay(rds@pVals,1)[1,1], 0.00522109, tolerance=1e-5)
   checkEqualsNumeric(assay(rds@pVals,1)[1,2], 0.01938319, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,4], 0.001032092, tolerance=1e-5)
}

# GRanges order
test_computeRoars_singleSamples_GRanges_order <- function() {
   gene_id <- c("B_POST", "A_PRE", "C_POST", "B_PRE", "A_POST", "C_PRE")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(c("-","+","+","-","+","+")),
      ranges = IRanges(
         start=c(20, 1, 52, 40, 10, 42),
         width=c(20, 10, 10, 5, 5, 10)),
      DataFrame(gene_id)
   )
   
   a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), cigar = "5M", strand = strand("+"))
   a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), cigar = "3M", strand = strand("+"))
   a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), cigar = "3M", strand = strand("+"))
   b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "1M", strand = strand("-"))
   b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), cigar = "3M", strand = strand("-"))
   # The next one is an overlapping read only if strandness is not considered. 
   # Otherwise it will be counted for pre_C. I will add a check on this in this test, even if it's not correct.
   overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), cigar = "5M", strand = strand("+"))
   c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), cigar = "50M", strand = strand("+"))
   
   
   treatmentAlign <- list(c(rep(a_pre, 2),rep(a_post, 3), a_pre_post, rep(b_pre,5), b_post, overlapbc))
   controlAlign <- list(c(a_post, rep(a_pre, 4), a_pre_post, rep(b_post,5), b_pre, c_post))
   
   rds <- RoarDataset(treatmentAlign, controlAlign, features)
   # Removed FALSE to test about default stranded value (FALSE).
   rds <- countPrePost(rds)
   rds <- computeRoars(rds)
   #assay(rds,2) <- as.matrix(data.frame(treatment_pre=mMtreatment, treatment_post=mMcontrol, control_pre=roar, control_post=pVal))
   checkEqualsNumeric(assay(rds,2)[1,1], -0.66538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,2], 1.21538461538462, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[1,3], -0.54746835443038, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[1,4]))
   checkEqualsNumeric(assay(rds,2)[2,1], 20.6923076923077, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,2], 0.04307692307692, tolerance=1e-5)
   checkEqualsNumeric(assay(rds,2)[2,3], 480.357142857144, tolerance=1e-5)
   checkTrue(is.na(assay(rds,2)[3,1]))
   checkEqualsNumeric(assay(rds,2)[3,2], -1)
   checkTrue(is.na(assay(rds,2)[3,3]))      
}
