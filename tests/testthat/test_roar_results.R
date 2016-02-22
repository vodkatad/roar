test_that("test_fpkmResults_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = c(1000, 900, 100, 2)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[1, 1] <- 5000
    assay(rds, 1)[1, 2] <- NA
    assay(rds, 1)[1, 3] <- 300
    assay(rds, 1)[1, 4] <- NA
    assay(rds, 1)[2, 1] <- 1000
    assay(rds, 1)[2, 2] <- NA
    assay(rds, 1)[2, 3] <- 10
    assay(rds, 1)[2, 4] <- NA
    dat <- fpkmResults(rds)
    expect_equal(833333, dat[1, "treatmentValue"], tolerance = .5, check.attributes=F)
    expect_equal(967742, dat[1, "controlValue"], tolerance = .5, check.attributes=F)
    expect_equal(1666667, dat[2, "treatmentValue"], tolerance = .5, check.attributes=F)
    expect_equal(322581, dat[2, "controlValue"], tolerance = .5, check.attributes=F)
})

test_that("test_fpkmResults_singleSamples_length", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    length <- c(1, 1, 50, 20)
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = c(1000, 900, 100, 2)), DataFrame(gene_id, 
            length))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[1, 1] <- 5000
    assay(rds, 1)[1, 2] <- NA
    assay(rds, 1)[1, 3] <- 300
    assay(rds, 1)[1, 4] <- NA
    assay(rds, 1)[2, 1] <- 1000
    assay(rds, 1)[2, 2] <- NA
    assay(rds, 1)[2, 3] <- 10
    assay(rds, 1)[2, 4] <- NA
    dat <- fpkmResults(rds)
    expect_equal(833333333, dat[1, "treatmentValue"], tolerance = .5, check.attributes=F)
    expect_equal(967741935, dat[1, "controlValue"], tolerance = .5, check.attributes=F)
    expect_equal(3333333, dat[2, "treatmentValue"], tolerance = .5, check.attributes=F)
    expect_equal(645161, dat[2, "controlValue"], tolerance = .5, check.attributes=F)
})

test_that("test_pvalueCorrectFilter_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "C_PRE", 
        "C_POST", "D_PRE", "D_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2, 3, 4, 5, 6), width = rep(1, length(gene_id))), 
        DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1, paired = FALSE)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 4, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[, "treatment_pre"] <- c(1, 1e+09, 1e+09, 20)
    assay(rds, 1)[, "control_pre"] <- c(20, 1e+09, 1e+09, 20)
    assay(rds, 2)[, "treatment_pre"] <- c(1, 1, 1, 1)
    assay(rds, 2)[, "treatment_post"] <- c(1, 1, 1, 1)
    assay(rds, 2)[, "control_pre"] <- c(1, 1, 1, 1)
    assay(rds, 2)[, "control_post"] <- c(2e-05, 0.001, 0.001, 
        0.01)
    dat <- pvalueCorrectFilter(rds, 1, 0.02, "bonferroni")
    dat_wanted <- data.frame(row.names = c("B", "C"), mM_treatment = c(1, 
        1), mM_control = c(1, 1), roar = c(1, 1), pval = c(0.003, 
        0.003), treatmentValue = c(5e+08, 5e+08), controlValue = c(5e+08, 
        5e+08))
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_pvalueFilter_mulSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = rep(1, length(gene_id))), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1, paired = FALSE)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 2
    length(rds@controlBams) <- 2
    assay(rds, 1)[, "treatment_pre"] <- c(1, 1e+09)
    assay(rds, 1)[, "control_pre"] <- c(20, 1e+09)
    assay(rds, 2)[, "treatment_pre"] <- c(1, 1)
    assay(rds, 2)[, "treatment_post"] <- c(2, 1)
    assay(rds, 2)[, "control_pre"] <- c(10, 1)
    assay(rds, 2)[, "control_post"] <- c(0.01, 0.1)
    rds@pVals <- SummarizedExperiment(assays = matrix(nrow = 2, 
        ncol = 4), rowRanges = preCoords, colData = DataFrame(row.names = c("pvalue_1_1", 
        "pvalue_1_2", "pvalue_2_1", "pvalue_2_2")))
    assay(rds@pVals, 1)[1, ] <- seq(1, 4)
    assay(rds@pVals, 1)[2, ] <- c(1, 1, 1, 5)
    dat <- pvalueFilter(rds, 0.01, 3)
    dat_wanted <- data.frame(row.names = c("A", "B"), mM_treatment = c(1, 
        1), mM_control = c(2, 1), roar = c(10, 1), pval = c(0.01, 
        0.1), pvalue_1_1 = c(1, 1), pvalue_1_2 = c(2, 1), pvalue_2_1 = c(3, 
        1), pvalue_2_2 = c(4, 5), treatmentValue = c(1, 1e+09), 
        controlValue = c(20, 1e+09), nUnderCutoff = c(2, 3))
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_pvalueFilter_mulSamples_paired", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = rep(1, length(gene_id))), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1, paired = TRUE)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 2
    length(rds@controlBams) <- 2
    assay(rds, 1)[, "treatment_pre"] <- c(1, 1e+09)
    assay(rds, 1)[, "control_pre"] <- c(20, 1e+09)
    assay(rds, 2)[, "treatment_pre"] <- c(1, 1)
    assay(rds, 2)[, "treatment_post"] <- c(2, 1)
    assay(rds, 2)[, "control_pre"] <- c(10, 1)
    assay(rds, 2)[, "control_post"] <- c(0.01, 0.1)
    rds@pVals <- SummarizedExperiment(assays = matrix(nrow = 2, 
        ncol = 2), rowRanges = preCoords, colData = DataFrame(row.names = c("pvalue_1_2", 
        "pvalue_2_1")))
    assay(rds@pVals, 1)[1, ] <- c(1, 20)
    assay(rds@pVals, 1)[2, ] <- c(5, 5)
    dat <- pvalueFilter(rds, 0.01, 0.05)
    dat_wanted <- data.frame(row.names = "A", mM_treatment = 1, 
        mM_control = 2, roar = 10, pval = 0.01, pvalue_1_2 = 1, 
        pvalue_2_1 = 20, treatmentValue = 1, controlValue = 20)
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_pvalueFilter_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "C_PRE", 
        "C_POST", "D_PRE", "D_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2, 3, 4, 5, 6), width = rep(1, length(gene_id))), 
        DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1, paired = FALSE)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 4, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[, "treatment_pre"] <- c(1, 1e+09, 1e+09, 20)
    assay(rds, 1)[, "control_pre"] <- c(20, 1e+09, 1e+09, 20)
    assay(rds, 2)[, "treatment_pre"] <- c(1, 1, -1, 2)
    assay(rds, 2)[, "treatment_post"] <- c(2, 1, 1, 1)
    assay(rds, 2)[, "control_pre"] <- c(10, 1, NA, 1)
    assay(rds, 2)[, "control_post"] <- c(0.01, 0.1, 0.1, 0.001)
    dat <- pvalueFilter(rds, 0.01, 0.05)
    dat_wanted <- data.frame(row.names = c("A", "D"), mM_treatment = c(1, 
        2), mM_control = c(2, 1), roar = c(10, 1), pval = c(0.01, 
        0.001), treatmentValue = c(0.5, 10), controlValue = c(10, 
        10))
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_standardFilter_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "C_PRE", 
        "C_POST", "D_PRE", "D_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2, 3, 4, 5, 6), width = rep(1, length(gene_id))), 
        DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1, paired = FALSE)
    preElems <- grep("_PRE$", mcols(rds@prePostCoords)$gene_id)
    postElems <- grep("_POST$", mcols(rds@prePostCoords)$gene_id)
    preCoords <- rds@prePostCoords[preElems, ]
    se <- SummarizedExperiment(assays = rep(list(matrix(nrow = 4, 
        ncol = 4)), 2), rowRanges = preCoords, colData = DataFrame(row.names = c("treatment_pre", 
        "treatment_post", "control_pre", "control_post")))
    rowRanges(rds) <- rowRanges(se)
    colData(rds) <- colData(se)
    assays(rds) <- assays(se)
    names(assays(rds)) <- "counts"
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 1)[, "treatment_pre"] <- c(1, 1e+09, 1e+09, 20)
    assay(rds, 1)[, "control_pre"] <- c(20, 1e+09, 1e+09, 20)
    assay(rds, 2)[, "treatment_pre"] <- c(1, -1, 1, 1)
    assay(rds, 2)[, "treatment_post"] <- c(1, 1, 1, 1)
    assay(rds, 2)[, "control_pre"] <- c(1, 1, NA, 1)
    assay(rds, 2)[, "control_post"] <- c(0.1, 0.1, 0.1, 0.1)
    dat <- standardFilter(rds, 1)
    dat_wanted <- data.frame(row.names = "D", mM_treatment = 1, 
        mM_control = 1, roar = 1, pval = 0.1, treatmentValue = 10, 
        controlValue = 10)
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_totalResults_mulSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = c(1000, 900, 1, 2)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 2
    length(rds@controlBams) <- 3
    assay(rds, 2)[1, 1] <- 1
    assay(rds, 2)[1, 2] <- 2
    assay(rds, 2)[1, 3] <- 3
    assay(rds, 2)[1, 4] <- 4
    assay(rds, 2)[2, 1] <- 5
    assay(rds, 2)[2, 2] <- 6
    assay(rds, 2)[2, 3] <- 7
    assay(rds, 2)[2, 4] <- 8
    rds@pVals <- SummarizedExperiment(assays = matrix(nrow = 2, 
        ncol = 6), rowRanges = preCoords, colData = DataFrame(row.names = c("pvalue_1_1", 
        "pvalue_1_2", "pvalue_1_3", "pvalue_2_1", "pvalue_2_2", 
        "pvalue_2_3")))
    assay(rds@pVals, 1)[1, ] <- seq(1, 6)
    assay(rds@pVals, 1)[2, ] <- seq(10, 15)
    dat <- totalResults(rds)
    dat_wanted <- data.frame(row.names = c("A", "B"), mM_treatment = c(1, 
        5), mM_control = c(2, 6), roar = c(3, 7), pval = c(4, 
        8), pvalue_1_1 = c(1, 10), pvalue_1_2 = c(2, 11), pvalue_1_3 = c(3, 
        12), pvalue_2_1 = c(4, 13), pvalue_2_2 = c(5, 14), pvalue_2_3 = c(6, 
        15))
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

test_that("test_totalResults_singleSamples", {
    gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
    features <- GRanges(seqnames = Rle(rep("chr1", length(gene_id))), 
        strand = strand(rep("+", length(gene_id))), ranges = IRanges(start = c(1000, 
            2000, 1, 2), width = c(1000, 900, 1, 2)), DataFrame(gene_id))
    se <- getPreCoordsSE(features)
    rds <- new("RoarDataset", se, treatmentBams = list(), controlBams = list(), 
        prePostCoords = features, step = 3, cores = 1)
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
    rds@postCoords <- rds@prePostCoords[postElems, ]
    length(rds@treatmentBams) <- 1
    length(rds@controlBams) <- 1
    assay(rds, 2)[1, 1] <- 1
    assay(rds, 2)[1, 2] <- 2
    assay(rds, 2)[1, 3] <- 3
    assay(rds, 2)[1, 4] <- 4
    assay(rds, 2)[2, 1] <- 5
    assay(rds, 2)[2, 2] <- 6
    assay(rds, 2)[2, 3] <- 7
    assay(rds, 2)[2, 4] <- 8
    dat <- totalResults(rds)
    dat_wanted <- data.frame(row.names = c("A", "B"), mM_treatment = c(1, 
        5), mM_control = c(2, 6), roar = c(3, 7), pval = c(4, 
        8))
    expect_equal(dat_wanted, dat, tolerance = 1e-6, check.attributes=F)
})

