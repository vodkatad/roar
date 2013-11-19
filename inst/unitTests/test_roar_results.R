# To manually perform single tests: source("inst/unitTests/test_roar_results.R")
# To perform all tests: BiocGenerics:::testPackage("roar")

test_totalResults_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2),
         width=c(1000, 900, 1, 2)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 1
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   assay(rds,2)[1,1] <- 1
   assay(rds,2)[1,2] <- 2
   assay(rds,2)[1,3] <- 3
   assay(rds,2)[1,4] <- 4
   assay(rds,2)[2,1] <- 5
   assay(rds,2)[2,2] <- 6
   assay(rds,2)[2,3] <- 7
   assay(rds,2)[2,4] <- 8
   
   dat <- totalResults(rds)
   dat_wanted <- data.frame(row.names=c("A","B"), 
                           mM_treatment=c(1,5), 
                           mM_control=c(2,6),
                           roar=c(3,7),
                           pval=c(4,8))
   checkEquals(dat, dat_wanted)
}

test_totalResults_mulSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2),
         width=c(1000, 900, 1, 2)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 2
   length(rds@controlBams)  <- 3
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   assay(rds,2)[1,1] <- 1
   assay(rds,2)[1,2] <- 2
   assay(rds,2)[1,3] <- 3
   assay(rds,2)[1,4] <- 4
   assay(rds,2)[2,1] <- 5
   assay(rds,2)[2,2] <- 6
   assay(rds,2)[2,3] <- 7
   assay(rds,2)[2,4] <- 8
   rds@pVals <- SummarizedExperiment(assays = matrix(nrow=2, ncol=6),
                                          rowData=preCoords, 
                                          colData=DataFrame(row.names=c("pvalue_1_1","pvalue_1_2","pvalue_1_3",
                                                                        "pvalue_2_1", "pvalue_2_2", "pvalue_2_3"))
   )
   assay(rds@pVals,1)[1,] <- seq(1,6)
   assay(rds@pVals,1)[2,] <- seq(10,15)
   
   dat <- totalResults(rds)
   dat_wanted <- data.frame(row.names=c("A","B"), 
                           mM_treatment=c(1,5), 
                           mM_control=c(2,6),
                           roar=c(3,7),
                           pval=c(4,8), 
                           pvalue_1_1 = c(1,10),
                           pvalue_1_2 = c(2,11),
                           pvalue_1_3 = c(3,12),
                           pvalue_2_1 = c(4,13),
                           pvalue_2_2 = c(5,14),
                           pvalue_2_3 = c(6,15))
   checkEquals(dat, dat_wanted)
}

test_fpkmResults_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2),
         width=c(1000, 900, 100, 2)),
      DataFrame(gene_id)
   )
   # prelen <- 1000, 100
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 1
   assay(rds,1)[1,1] <- 5000 # treatment_pre A
   assay(rds,1)[1,2] <- NA
   assay(rds,1)[1,3] <- 300 # control_pre A
   assay(rds,1)[1,4] <- NA
   assay(rds,1)[2,1] <- 1000 # treatment_pre B
   assay(rds,1)[2,2] <- NA
   assay(rds,1)[2,3] <- 10 # control_pre B
   assay(rds,1)[2,4] <- NA
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   
   dat <- fpkmResults(rds)

   checkEqualsNumeric(dat[1,"treatmentValue"], 8.3333e5, tolerance=1e-5)
   checkEqualsNumeric(dat[1,"controlValue"], 9.67741e5, tolerance=1e-5)
   checkEqualsNumeric(dat[2,"treatmentValue"], 1.666667e6, tolerance=1e-5)
   checkEqualsNumeric(dat[2,"controlValue"], 3.2258e5, tolerance=1e-5)
}

test_standardFilter_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "C_PRE", "C_POST", "D_PRE", "D_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2, 3, 4, 5, 6),
         width=rep(1, length(gene_id))),
      DataFrame(gene_id)
   )
   # prelen A-1, B-1, C-1, D-1
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=4, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 1
   # Here we need to test filter on: FPKM (it will be A). negative m/M (B). NA roar (C). Only D will survive.
   # For A and D we will need FPKM values, <= 1 for A and > 1 for D.
   # We set pre counts for all genes:
   assay(rds,1)[,"treatment_pre"] <- c(1,1e9,1e9,20)
   assay(rds,1)[,"control_pre"] <- c(20,1e9,1e9,20)
   # D treatment FPKM 9.9999
   
   # We set a negative m/M for B, a NA roar for C. Other values are ok.
   assay(rds,2) <- matrix(ncol=4, nrow=4)
   assay(rds,2)[, "treatment_pre"] <- c(1,-1,1,1)
   assay(rds,2)[, "treatment_post"] <- c(1,1,1,1)
   assay(rds,2)[, "control_pre"] <- c(1,1,NA,1)
   assay(rds,2)[, "control_post"] <- c(0.1,0.1,0.1,0.1)
      
   dat <- standardFilter(rds, 1)
   dat_wanted <- data.frame(row.names="D", 
                           mM_treatment=1, 
                           mM_control=1,
                           roar=1,
                           pval=0.1,
                           treatmentValue=10,
                           controlValue=10)
                           #bonferroniPval=0.1)
   checkEquals(dat, dat_wanted, tolerance=1e-5)
   # XXX TODO check whi 9.9999 is rounded to 10
}
   
test_pvalueFilter_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "C_PRE", "C_POST", "D_PRE", "D_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2, 3, 4, 5, 6),
         width=rep(1, length(gene_id))),
      DataFrame(gene_id)
   )
   # prelen A-1, B-1, C-1, D-1
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=4, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 1
   length(rds@controlBams)  <- 1
   # Here we need to test filter on pvalues, therefore we will need three genes (A, B, C)
   # passing the FPKM/mM/roar (A/B/D) filter and then two (A/D) of them passing the pvalue filter.
   # We set pre counts for all genes:
   assay(rds,1)[,"treatment_pre"] <- c(1,1e9,1e9,20)
   assay(rds,1)[,"control_pre"] <- c(20,1e9,1e9,20)
   # D treatment FPKM 9.9999
   
   # We set a negative m/M for C (and a NA roar). Other values are ok.
   assay(rds,2) <- matrix(ncol=4, nrow=4)
   assay(rds,2)[, "treatment_pre"] <- c(1,1,-1,2)
   assay(rds,2)[, "treatment_post"] <- c(2,1,1,1)
   assay(rds,2)[, "control_pre"] <- c(10,1,NA,1)
   assay(rds,2)[, "control_post"] <- c(0.01,0.1,0.1,0.001)
   # The pvalues will be bonferroni corrected with 3 genes. Then only 2 will survive the pvalue filter.
   
   dat <- pvalueFilter(rds, 0.01, 0.05)
   dat_wanted <- data.frame(row.names=c("A","D"), 
                           mM_treatment=c(1,2), 
                           mM_control=c(2,1),
                           roar=c(10,1),
                           pval=c(0.01, 0.001),
                           treatmentValue=c(0.5, 10),
                           controlValue=c(10, 10))
                           #bonferroniPval=c(0.03, 0.003))
   checkEquals(dat, dat_wanted, tolerance=1e-5)
   # XXX TODO check whi 9.9999 is rounded to 10
}

test_pvalueFilter_mulSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(rep("chr1", length(gene_id))),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 1, 2),
         width=rep(1, length(gene_id))),
      DataFrame(gene_id)
   )
   # prelen A-1, B-1
   rds <- new("RoarDataset", treatmentBams=list(), controlBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@treatmentBams)  <- 2
   length(rds@controlBams)  <- 2
   assay(rds,1)[,"treatment_pre"] <- c(1,1e9)
   assay(rds,1)[,"control_pre"] <- c(20,1e9)
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   assay(rds,2)[, "treatment_pre"] <- c(1,1)
   assay(rds,2)[, "treatment_post"] <- c(2,1)
   assay(rds,2)[, "control_pre"] <- c(10,1)
   assay(rds,2)[, "control_post"] <- c(0.01,0.1)
   # To test multiple samples we keep all the values equal to the previous test (with only two genes), 
   # therefore we will have a working object. We just need to add the rds@pVals slot with the needed data.
   rds@pVals <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                                     rowData=preCoords, 
                                     colData=DataFrame(row.names=c("pvalue_1_1","pvalue_1_2","pvalue_2_1",
                                                                   "pvalue_2_2"))
   )
   assay(rds@pVals,1)[1,] <- seq(1,4)
   assay(rds@pVals,1)[2,] <- c(1,1,1,5)
   
   
   dat <- pvalueFilter(rds, 0.01, 3)
   dat_wanted <- data.frame(row.names=c("A","B"), 
                           mM_treatment=c(1,1), 
                           mM_control=c(2,1),
                           roar=c(10,1),
                           pval=c(0.01, 0.1),
                           pvalue_1_1 = c(1,1),
                           pvalue_1_2 = c(2,1),
                           pvalue_2_1 = c(3,1),
                           pvalue_2_2 = c(4,5),
                           treatmentValue=c(1, 1e9),
                           controlValue=c(20, 1e9),
                           #bonferroniPval=c(0.02, 0.2)
                           nUnderCutoff=c(2,3))
   checkEquals(dat, dat_wanted, tolerance=1e-5)
   # XXX TODO check whi 9.9999 is rounded to 10
}
