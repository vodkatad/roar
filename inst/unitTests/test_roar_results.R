# To manually perform single tests: source("inst/unitTests/test_roar.R")

test_totalResults_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000),
         width=c(1000, 900)),
      DataFrame(gene_id)
   )
   rds <- new("RoarDataset", rightBams=list(), leftBams=list(), 
              prePostCoords=features, step = 2, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=1, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
   )
   
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   # We need to set these lengths to choose the right branch of the if in computePvals (about number of samples).
   length(rds@rightBams)  <- 3
   length(rds@leftBams)  <- 2
   rds@countsRight <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                           rowData=preCoords, 
                                           colData=DataFrame(row.names=c("pre","post"))
   )
   rds@countsLeft <- SummarizedExperiment(assays = matrix(nrow=1, ncol=2),
                                          rowData=preCoords, 
                                          colData=DataFrame(row.names=c("pre","post"))
   )
   assay(rds@countsRight,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsRight,3) <- matrix(nrow=1,ncol=2)
   assay(rds@countsLeft,2) <- matrix(nrow=1,ncol=2)
   assay(rds@countsRight,1)[1,] <- c(10,10)
   assay(rds@countsRight,2)[1,] <- c(1,10)
   assay(rds@countsRight,3)[1,] <- c(20,10)
   assay(rds@countsLeft,1)[1,] <- c(10,5)
   assay(rds@countsLeft,2)[1,] <- c(10,20)
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
