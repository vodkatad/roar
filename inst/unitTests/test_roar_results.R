# To manually perform single tests: source("inst/unitTests/test_roar_results.R")

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
   rds <- new("RoarDataset", rightBams=list(), leftBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@rightBams)  <- 1
   length(rds@leftBams)  <- 1
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   assay(rds,2)[1,1] <- 1
   assay(rds,2)[1,2] <- 2
   assay(rds,2)[1,3] <- 3
   assay(rds,2)[1,4] <- 4
   assay(rds,2)[2,1] <- 5
   assay(rds,2)[2,2] <- 6
   assay(rds,2)[2,3] <- 7
   assay(rds,2)[2,4] <- 8
   
   df <- totalResults(rds)
   df_wanted <- data.frame(row.names=c("A","B"), 
                           mM_right=c(1,5), 
                           mM_left=c(2,6),
                           roar=c(3,7),
                           pval=c(4,8))
   checkEquals(df, df_wanted)
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
   rds <- new("RoarDataset", rightBams=list(), leftBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@rightBams)  <- 2
   length(rds@leftBams)  <- 3
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
   
   df <- totalResults(rds)
   df_wanted <- data.frame(row.names=c("A","B"), 
                           mM_right=c(1,5), 
                           mM_left=c(2,6),
                           roar=c(3,7),
                           pval=c(4,8), 
                           pvalue_1_1 = c(1,10),
                           pvalue_1_2 = c(2,11),
                           pvalue_1_3 = c(3,12),
                           pvalue_2_1 = c(4,13),
                           pvalue_2_2 = c(5,14),
                           pvalue_2_3 = c(6,15))
   checkEquals(df, df_wanted)
}


test_filteringInfoResults_singleSamples <- function() {
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
   rds <- new("RoarDataset", rightBams=list(), leftBams=list(), 
              prePostCoords=features, step = 3, cores=1)
   preElems <- grep("_PRE$", elementMetadata(rds@prePostCoords)$gene_id)
   postElems <- grep("_POST$", elementMetadata(rds@prePostCoords)$gene_id)
   preCoords <- rds@prePostCoords[preElems,]
   se <- SummarizedExperiment(assays = matrix(nrow=2, ncol=4),
                              rowData=preCoords, 
                              colData=DataFrame(row.names=c("right_pre","right_post","left_pre", "left_post"))
   )
   rowData(rds) <- rowData(se)
   colData(rds) <- colData(se)
   assays(rds) <- assays(se)
   names(assays(rds)) <- "counts"
   rds@postCoords <- rds@prePostCoords[postElems,]
   length(rds@rightBams)  <- 1
   length(rds@leftBams)  <- 1
   assay(rds,1)[1,1] <- 5000 # right_pre A
   assay(rds,1)[1,2] <- NA
   assay(rds,1)[1,3] <- 300 # left_pre A
   assay(rds,1)[1,4] <- NA
   assay(rds,1)[2,1] <- 1000 # right_pre B
   assay(rds,1)[2,2] <- NA
   assay(rds,1)[2,3] <- 10 # left_pre B
   assay(rds,1)[2,4] <- NA
   assay(rds,2) <- matrix(ncol=4, nrow=2)
   
   df <- filteringInfoResults(rds)

   checkEqualsNumeric(df[1,"rightValue"], 8.3333e5, tolerance=1e-5)
   checkEqualsNumeric(df[1,"leftValue"], 9.67741e5, tolerance=1e-5)
   checkEqualsNumeric(df[2,"rightValue"], 1.666667e6, tolerance=1e-5)
   checkEqualsNumeric(df[2,"leftValue"], 3.2258e5, tolerance=1e-5)
}