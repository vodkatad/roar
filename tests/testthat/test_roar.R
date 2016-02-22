# Bad! What should I do? Export this function is not an option.
getPreCoordsSE <- function(gtfGRanges) {
   # Now we need to keep means and totals of counts over PRE/POST for the two lists.
   # In the simpler case with a single alignment for both conditions we just keep the counts.
   preElems <- grep("_PRE$", mcols(gtfGRanges)$gene_id)
   preCoords <- gtfGRanges[preElems,]
   se <- SummarizedExperiment(rowRanges=preCoords, colData=DataFrame(row.names=c("treatment_pre","treatment_post","control_pre", "control_post")))
   return(se)
}

test_that("test_computePairedPvals", {
    gene_id <- c("A_PRE", "A_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1")), strand = strand(rep("+", 
        length(gene_id))), ranges = IRanges(start = c(1000, 2000), 
        width = c(1000, 900)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 2, cores = 1)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    length(rds@treatmentBams) <- 3
    length(rds@controlBams) <- 2
    rds@countsTreatment <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 2)), 3), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    rds@countsControl <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 2)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    assay(rds@countsTreatment, 1)[1, ] <- c(10, 10)
    assay(rds@countsTreatment, 2)[1, ] <- c(1, 10)
    assay(rds@countsTreatment, 3)[1, ] <- c(20, 10)
    assay(rds@countsControl, 1)[1, ] <- c(10, 5)
    assay(rds@countsControl, 2)[1, ] <- c(10, 20)
    rds <- computePairedPvals(rds, c(2, 3), c(1, 2))
    expect_equal(0.00522109, assay(rds@pVals, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.01938319, assay(rds@pVals, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0.001032092, assay(rds, 2)[1, 4], tolerance=1e-6, check.attributes=F)
})

test_that("test_computePvals_multipleSamples", {
    gene_id <- c("A_PRE", "A_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1")), strand = strand(rep("+", 
        length(gene_id))), ranges = IRanges(start = c(1000, 2000), 
        width = c(1000, 900)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 2, cores = 1)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    length(rds@treatmentBams) <- 3
    length(rds@controlBams) <- 2
    rds@countsTreatment <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 2)), 3), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    rds@countsControl <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 2)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    assay(rds@countsTreatment, 1)[1, ] <- c(10, 10)
    assay(rds@countsTreatment, 2)[1, ] <- c(1, 10)
    assay(rds@countsTreatment, 3)[1, ] <- c(20, 10)
    assay(rds@countsControl, 1)[1, ] <- c(10, 5)
    assay(rds@countsControl, 2)[1, ] <- c(10, 20)
    rds <- computePvals(rds)
    expect_equal(0.491599, assay(rds@pVals, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.257549, assay(rds@pVals, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0.00522109, assay(rds@pVals, 1)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(0.23264154, assay(rds@pVals, 1)[1, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds@pVals, 1)[1, 5], tolerance=1e-6, check.attributes=F)
    expect_equal(0.01938319, assay(rds@pVals, 1)[1, 6], tolerance=1e-6, check.attributes=F)
})

test_that("test_computePvals_single", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2", 
        "chr2")), strand = strand(rep("+", length(gene_id))), 
        ranges = IRanges(start = c(1000, 2000, 3000, 3600), width = c(1000, 
            900, 600, 300)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 2, cores = 1)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 2, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[1, 1] <- 10
    assay(rds, 1)[1, 2] <- 100
    assay(rds, 1)[1, 3] <- 20
    assay(rds, 1)[1, 4] <- 20
    assay(rds, 1)[2, 1] <- 10
    assay(rds, 1)[2, 2] <- 20
    assay(rds, 1)[2, 3] <- 20
    assay(rds, 1)[2, 4] <- 20
    rds <- computePvals(rds)
    expect_equal(2.212406e-07, assay(rds, 2)[1, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(0.22337243, assay(rds, 2)[2, 4], tolerance=1e-6, check.attributes=F)
})

test_that("test_computePvals_singlevsMul", {
    gene_id <- c("A_PRE", "A_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1")), strand = strand(rep("+", 
        length(gene_id))), ranges = IRanges(start = c(1000, 2000), 
        width = c(1000, 900)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 2, cores = 1)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 3
    rds@countsTreatment <- SummarizedExperiment(assays = matrix(nrow = 1, 
        ncol = 2), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    rds@countsControl <- SummarizedExperiment(assays = rep(list(matrix(nrow = 1, 
        ncol = 2)), 3), rowRanges = preCoords, colData = DataFrame(row.names = c("pre", 
        "post")))
    assay(rds@countsTreatment, 1)[1, ] <- c(10, 10)
    assay(rds@countsControl, 1)[1, ] <- c(10, 10)
    assay(rds@countsControl, 2)[1, ] <- c(10, 5)
    assay(rds@countsControl, 3)[1, ] <- c(10, 20)
    rds <- computePvals(rds)
    expect_equal(1, assay(rds@pVals, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.491599, assay(rds@pVals, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0.257549, assay(rds@pVals, 1)[1, 3], tolerance=1e-6, check.attributes=F)
})

test_that("test_computeRoars_multipleSamples", {
    gene_id <- c("A_PRE", "A_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1")), strand = strand(c("+", 
        "+")), ranges = IRanges(start = c(1, 10), width = c(10, 
        5)), DataFrame(gene_id))
    a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
        cigar = "5M", strand = strand("+"))
    a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), 
        cigar = "3M", strand = strand("+"))
    treatmentAlign <- list(rep(a_pre, 2), rep(a_post, 3))
    controlAlign <- list(a_post, rep(a_pre, 4), a_post)
    rds <- RoarDataset(treatmentAlign, controlAlign, features)
    rds <- countPrePost(rds, FALSE)
    rds <- computeRoars(rds)
    expect_equal(-0.48, assay(rds, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.666666666, assay(rds, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(-0.72, assay(rds, 2)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[1, 4]))
})

test_that("test_computeRoars_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_POST", "B_PRE", "C_PRE", 
        "C_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(c("+", "+", "-", "-", "+", "+")), ranges = IRanges(start = c(1, 
            10, 20, 40, 42, 52), width = c(10, 5, 20, 5, 10, 
            10)), DataFrame(gene_id))
    a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
        cigar = "5M", strand = strand("+"))
    a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), 
        cigar = "3M", strand = strand("+"))
    a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), 
        cigar = "3M", strand = strand("+"))
    b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "1M", strand = strand("-"))
    b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), 
        cigar = "3M", strand = strand("-"))
    overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "5M", strand = strand("+"))
    c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), 
        cigar = "50M", strand = strand("+"))
    treatmentAlign <- list(c(rep(a_pre, 2), rep(a_post, 3), a_pre_post, 
        rep(b_pre, 5), b_post, overlapbc))
    controlAlign <- list(c(a_post, rep(a_pre, 4), a_pre_post, 
        rep(b_post, 5), b_pre, c_post))
    rds <- RoarDataset(treatmentAlign, controlAlign, features)
    rds <- countPrePost(rds, FALSE)
    rds <- computeRoars(rds)
    expect_equal(-0.66538461538462, assay(rds, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(1.21538461538462, assay(rds, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(-0.54746835443038, assay(rds, 2)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[1, 4]))
    expect_equal(20.6923076923077, assay(rds, 2)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.04307692307692, assay(rds, 2)[2, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(480.357142857144, assay(rds, 2)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 1]))
    expect_equal(-1, assay(rds, 2)[3, 2], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 3]))
    rds2 <- RoarDataset(treatmentAlign, controlAlign, features)
    rds2 <- countPrePost(rds2, TRUE)
    expect_equal(2, assay(rds2, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(4, assay(rds2, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(4, assay(rds2, 1)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(2, assay(rds2, 1)[1, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(5, assay(rds2, 1)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds2, 1)[2, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds2, 1)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(5, assay(rds2, 1)[2, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds2, 1)[3, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds2, 1)[3, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds2, 1)[3, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds2, 1)[3, 4], tolerance=1e-6, check.attributes=F)
})

test_that("test_computeRoars_singleSamples_GRanges_order", {
    gene_id <- c("B_POST", "A_PRE", "C_POST", "B_PRE", "A_POST", 
        "C_PRE")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(c("-", "+", "+", "-", "+", "+")), ranges = IRanges(start = c(20, 
            1, 52, 40, 10, 42), width = c(20, 10, 10, 5, 5, 10)), 
        DataFrame(gene_id))
    a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
        cigar = "5M", strand = strand("+"))
    a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), 
        cigar = "3M", strand = strand("+"))
    a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), 
        cigar = "3M", strand = strand("+"))
    b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "1M", strand = strand("-"))
    b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), 
        cigar = "3M", strand = strand("-"))
    overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "5M", strand = strand("+"))
    c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), 
        cigar = "50M", strand = strand("+"))
    treatmentAlign <- list(c(rep(a_pre, 2), rep(a_post, 3), a_pre_post, 
        rep(b_pre, 5), b_post, overlapbc))
    controlAlign <- list(c(a_post, rep(a_pre, 4), a_pre_post, 
        rep(b_post, 5), b_pre, c_post))
    rds <- RoarDataset(treatmentAlign, controlAlign, features)
    rds <- countPrePost(rds)
    rds <- computeRoars(rds)
    expect_equal(-0.66538461538462, assay(rds, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(1.21538461538462, assay(rds, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(-0.54746835443038, assay(rds, 2)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[1, 4]))
    expect_equal(20.6923076923077, assay(rds, 2)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.04307692307692, assay(rds, 2)[2, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(480.357142857144, assay(rds, 2)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 1]))
    expect_equal(-1, assay(rds, 2)[3, 2], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 3]))
})

test_that("test_computeRoars_singleSamples_length", {
    gene_id <- c("A_PRE", "A_POST", "B_POST", "B_PRE", "C_PRE", 
        "C_POST")
    length <- c(1, 15, 20, 5, 30, 40)
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(c("+", "+", "-", "-", "+", "+")), ranges = IRanges(start = c(1, 
            10, 20, 40, 42, 52), width = c(10, 5, 20, 5, 10, 
            10)), DataFrame(gene_id, length))
    a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
        cigar = "5M", strand = strand("+"))
    a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), 
        cigar = "3M", strand = strand("+"))
    a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), 
        cigar = "3M", strand = strand("+"))
    b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "1M", strand = strand("-"))
    b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), 
        cigar = "3M", strand = strand("-"))
    overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "5M", strand = strand("+"))
    c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), 
        cigar = "50M", strand = strand("+"))
    treatmentAlign <- list(c(rep(a_pre, 2), rep(a_post, 3), a_pre_post, 
        rep(b_pre, 5), b_post, overlapbc))
    controlAlign <- list(c(a_post, rep(a_pre, 4), a_pre_post, 
        rep(b_post, 5), b_pre, c_post))
    rds <- RoarDataset(treatmentAlign, controlAlign, features)
    rds <- countPrePost(rds, FALSE)
    rds <- computeRoars(rds)
    expect_equal(7.34615384615385, assay(rds, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(41.1538461538462, assay(rds, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0.178504672897196, assay(rds, 2)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[1, 4]))
    expect_equal(20.6923076923077, assay(rds, 2)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.04307692307692, assay(rds, 2)[2, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(480.357142857144, assay(rds, 2)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 1]))
    expect_equal(-1, assay(rds, 2)[3, 2], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 3]))
})

test_that("test_computeRoars_singlevsMulSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_POST", "B_PRE", "C_PRE", 
        "C_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(c("+", "+", "-", "-", "+", "+")), ranges = IRanges(start = c(1, 
            10, 20, 40, 42, 52), width = c(10, 5, 20, 5, 10, 
            10)), DataFrame(gene_id))
    a_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
        cigar = "5M", strand = strand("+"))
    a_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(11), 
        cigar = "3M", strand = strand("+"))
    a_pre_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(9), 
        cigar = "3M", strand = strand("+"))
    b_pre <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "1M", strand = strand("-"))
    b_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(21), 
        cigar = "3M", strand = strand("-"))
    overlapbc <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(40), 
        cigar = "5M", strand = strand("+"))
    c_post <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(53), 
        cigar = "50M", strand = strand("+"))
    treatmentAlign <- list(c(rep(a_pre, 2), rep(a_post, 3), a_pre_post, 
        rep(b_pre, 5), b_post, overlapbc))
    controlAlign <- list(a_post, rep(a_pre, 4), a_pre_post, rep(b_post, 
        5), b_pre, c_post)
    rds <- RoarDataset(treatmentAlign, controlAlign, features)
    rds <- countPrePost(rds, FALSE)
    rds <- computeRoars(rds)
    expect_equal(-0.66538461538462, assay(rds, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(1.21538461538462, assay(rds, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(-0.54746835443038, assay(rds, 2)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[1, 4]))
    expect_equal(20.6923076923077, assay(rds, 2)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0.04307692307692, assay(rds, 2)[2, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(480.357142857144, assay(rds, 2)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 1]))
    expect_equal(-1, assay(rds, 2)[3, 2], tolerance=1e-6, check.attributes=F)
    expect_true(is.na(assay(rds, 2)[3, 3]))
})

test_that("test_countPrePost_mulSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", 
        "D_POST", "E_PRE", "E_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2", 
        "chr2", "chr2", "chr2", "chr1", "chr1")), strand = strand(rep("+", 
        length(gene_id))), ranges = IRanges(start = c(1000, 2000, 
        3000, 3600, 7000, 7500, 4000, 4000), width = c(500, 900, 
        500, 300, 600, 300, 500, 900)), DataFrame(gene_id))
    rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), 
        cigar = "300M", strand = strand("+"))
    rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), 
        cigar = "300M", strand = strand("+"))
    rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), 
        cigar = "300M", strand = strand("+"))
    rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), 
        cigar = "300M", strand = strand("+"))
    rds <- RoarDataset(list(rd1, rd2), list(rd3, rd4), features)
    rds <- countPrePost(rds, FALSE)
    expect_equal(1, assay(rds@countsTreatment, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds@countsTreatment, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds@countsControl, 2)[2, 1], tolerance=1e-6, check.attributes=F)
})

test_that("test_countPrePost_preferPOST", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2", 
        "chr2")), strand = strand(rep("+", length(gene_id))), 
        ranges = IRanges(start = c(1000, 2000, 3000, 3600), width = c(1000, 
            900, 600, 300)), DataFrame(gene_id))
    rd1 <- GAlignments("a", seqnames = factor("chr1", levels = c("chr1", 
        "chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
    rd2 <- GAlignments("a", seqnames = factor("chr1", levels = c("chr1", 
        "chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
    rd3 <- GAlignments("a", seqnames = factor("chr2", levels = c("chr1", 
        "chr2")), pos = as.integer(2800), cigar = "300M", strand = strand("+"))
    rd4 <- GAlignments("a", seqnames = factor("chr2", levels = c("chr1", 
        "chr2")), pos = as.integer(3500), cigar = "300M", strand = strand("+"))
    rd5 <- GAlignments("a", seqnames = factor("chr1", levels = c("chr1", 
        "chr2")), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
    rds <- RoarDataset(list(c(rd1, rd2)), list(c(rd3, rd4, rd5)), 
        features)
    rds <- countPrePost(rds, FALSE)
    expect_equal(0, assay(rds, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(2, assay(rds, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds, 1)[1, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds, 1)[1, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds, 1)[2, 3], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds, 1)[2, 4], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds, 1)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds, 1)[2, 2], tolerance=1e-6, check.attributes=F)
})

test_that("test_countPrePost_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", 
        "D_POST", "E_PRE", "E_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2", 
        "chr2", "chr2", "chr2", "chr1", "chr1")), strand = strand(rep("+", 
        length(gene_id))), ranges = IRanges(start = c(1000, 2000, 
        3000, 3600, 7000, 7500, 4000, 4000), width = c(500, 900, 
        500, 300, 600, 300, 500, 900)), DataFrame(gene_id))
    rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), 
        cigar = "300M", strand = strand("+"))
    rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), 
        cigar = "300M", strand = strand("+"))
    rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), 
        cigar = "300M", strand = strand("+"))
    rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), 
        cigar = "300M", strand = strand("+"))
    rds <- RoarDataset(list(c(rd1, rd2)), list(c(rd3, rd4)), 
        features)
    rds <- countPrePost(rds, FALSE)
    expect_equal(1, assay(rds, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(2, assay(rds, 1)[2, 3], tolerance=1e-6, check.attributes=F)
})

test_that("test_countPrePost_stranded", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(c("chr1", "chr1", "chr2", 
        "chr2")), strand = strand(rep("+", length(gene_id))), 
        ranges = IRanges(start = c(1000, 2000, 3000, 3600), width = c(1000, 
            900, 600, 300)), DataFrame(gene_id))
    rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), 
        cigar = "300M", strand = strand("+"))
    rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), 
        cigar = "300M", strand = strand("-"))
    rd3 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), 
        cigar = "300M", strand = strand("+"))
    rd4 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), 
        cigar = "300M", strand = strand("-"))
    rd5 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), 
        cigar = "300M", strand = strand("+"))
    rd6 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), 
        cigar = "300M", strand = strand("+"))
    rds <- RoarDataset(list(rd1, rd2), list(rd3, rd4, c(rd5, 
        rd6)), features)
    rds <- countPrePost(rds, TRUE)
    expect_equal(1, assay(rds@countsTreatment, 1)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds@countsTreatment, 1)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds@countsTreatment, 2)[1, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds@countsTreatment, 2)[1, 2], tolerance=1e-6, check.attributes=F)
    expect_equal(1, assay(rds@countsControl, 1)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(0, assay(rds@countsControl, 2)[2, 1], tolerance=1e-6, check.attributes=F)
    expect_equal(2, assay(rds@countsControl, 3)[2, 1], tolerance=1e-6, check.attributes=F)
})

