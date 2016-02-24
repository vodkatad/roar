test_that("test_getApaGenesFractionsPlusStrandBasic", {
    # Here we put two APA inside the first and last exon of a three exon gene.
    chr <- rep("chr1", 3)
    strand <- rep("+", 3)
    gene <- rep("A", 3)
    type <- rep("gene", 3)
    apa <- rep(NA, 3)
    geneGr <- GRanges(seqnames = chr, strand = strand, ranges = IRanges(start = c(1000, 
        2000, 3300), width = c(500, 900, 10)), DataFrame(gene, 
        apa, type))
    apas <- GRanges(seqnames = c(chr[1], chr[1]), strand = c(strand[1], 
        strand[1]), ranges = IRanges(start = c(1100, 3305), width = c(1, 
        1)), DataFrame(gene = rep(gene[1], 2), apa = paste0(c("apa1_", 
        "apa2_"), gene[1]), type = rep("apa", 2)))
    res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr[1], 
        strand[1], gene[1])
    # Two possible apa choices.
    expect_equal(2, length(res)) 
    # Four fragments.
    expect_equal(4, length(res[[1]]))
    # Two info for which fragments need to be considered for the two choices.
    expect_equal(2, length(res[[2]]))
    # For the first APA choice we sum for PRE only the first fragment counts.
    expect_equal(1, res[[2]][[1]]@PREstart)
    expect_equal(1, res[[2]][[1]]@PREend)
    # For the second one up to the third (for POST we get all the remaining fragments counts).
    expect_equal(3, res[[2]][[2]]@PREstart)
    expect_equal(3, res[[2]][[2]]@PREstart)
    expect_equal("apa2", res[[2]][[2]]@name)
    expect_equal("apa1", res[[2]][[1]]@name)
    wanted <- GRanges(seqnames = rep(chr[1], 4), strand = rep(strand[1], 
        4), ranges = IRanges(start = c(1000, 1101, 3300, 3306), 
        width = c(101, 2199, 6, 4)), DataFrame(length = c(101, 
        1299, 6, 4)))
    expect_identical(res[[1]], wanted)
})

test_that("test_getApaGenesFractionsPlusStrandIntronExon", {
    chr <- "chr1"
    strand <- "+"
    gene <- "A"
    type <- "gene"
    apa <- NA
    geneGr <- GRanges(seqnames = rep(chr, 3), strand = rep(strand, 
        3), ranges = IRanges(start = c(1000, 2000, 3300), width = c(500, 
        900, 10)), DataFrame(rep(gene, 3), rep(apa, 3), rep(type, 
        3)))
    apas <- GRanges(seqnames = rep(chr, 2), strand = rep(strand, 
        2), ranges = IRanges(start = c(2010, 2901), width = c(1, 
        1)), DataFrame(gene = rep(gene, 2), apa = paste0(c("apa1_", 
        "apa2_"), gene), type = rep("apa", 2)))
    res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr, 
        strand, gene)
    expect_equal(4, length(res[[1]]))
    expect_equal(2, length(res[[2]]))
    expect_equal(1, res[[2]][[1]]@PREstart)
    expect_equal(1, res[[2]][[1]]@PREend)
    expect_equal(1, res[[2]][[2]]@PREstart)
    expect_equal(2, res[[2]][[2]]@PREend)
    wanted <- GRanges(seqnames = rep(chr[1], 4), strand = rep(strand[1], 
        4), ranges = IRanges(start = c(2000, 2011, 2902, 3300), 
        width = c(11, 891, 398, 10)), DataFrame(length = c(11, 
        889, 0, 10)))
    expect_identical(res[[1]], wanted)
})

test_that("test_getApaGenesFractionsPlusStrandSingleIntron", 
    {
        chr <- "chr1"
        strand <- "+"
        gene <- "A"
        type <- "gene"
        apa <- NA
        geneGr <- GRanges(seqnames = rep(chr, 3), strand = rep(strand, 
            3), ranges = IRanges(start = c(1000, 2000, 3300), 
            width = c(500, 900, 10)), DataFrame(rep(gene, 3), 
            rep(apa, 3), rep(type, 3)))
        apas <- GRanges(seqnames = chr, strand = strand, ranges = IRanges(start = 1510, 
            width = 1), DataFrame(gene = gene, apa = paste0("apa1_", 
            gene), type = "apa"))
        res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr, 
            strand, gene)
        expect_equal(3, length(res[[1]]))
        expect_equal(1, length(res[[2]]))
        expect_equal(1, res[[2]][[1]]@PREstart)
        expect_equal(1, res[[2]][[1]]@PREend)
        expect_equal("apa1", res[[2]][[1]]@name)
        wanted <- GRanges(seqnames = rep(chr, 3), strand = rep(strand, 
            3), ranges = IRanges(start = c(1000, 1511, 3300), 
            width = c(511, 1789, 10)), DataFrame(length = c(500, 
            900, 10)))
        expect_identical(res[[1]], wanted)
    })

test_that("test_RoarDatasetMultipleAPA_error_mismatch_gene_apa", 
    {
        gene <- c("A", "B", NA)
        type <- c("gene", "gene", "apa")
        apa <- c(NA, NA, "apa1_A")
        features <- GRanges(seqnames = Rle(c("chr1", "chr2", 
            "chr1")), strand = strand(rep("+", length(gene))), 
            ranges = IRanges(start = c(1000, 2000, 1300), width = c(500, 
                900, 1)), DataFrame(gene, apa, type))
        rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), 
            cigar = "300M", strand = strand("+"))
        obs <- tryCatch(RoarDatasetMultipleAPA(list(c(rd1, rd1)), 
            list(c(rd1, rd1)), features), error = function(e) e)
        expect_equal("All the genes in the gtf should have at least one apa", 
            obs$message)
    })

test_that("test_badIntegrationTest", 
    {
       chr <- rep("chr1", 3)
       strand <- rep("+", 3)
       gene <- rep("A", 3)
       type <- rep("gene", 3)
       apa <- rep(NA, 3)
       geneGr <- GRanges(seqnames = chr, strand = strand, ranges = IRanges(start = c(1000, 
                                                                                     2000, 3300), width = c(500, 900, 10)), DataFrame(gene, apa, type))
       gene2 <- rep("B", 3)
       chr2 = rep("chr2", 3)
       geneGr2 <- GRanges(seqnames = chr2, strand = strand, ranges = IRanges(start = c(1000, 
                                                                                       2000, 3300), width = c(500, 900, 10)), DataFrame(gene=gene2, apa, type))
       apas1 <- GRanges(seqnames = c(chr[1], chr[1]), strand = c(strand[1], strand[1]), 
                        ranges = IRanges(start = c(1100, 3305), width = c(1, 1)), 
                        DataFrame(gene = rep(gene[1], 2), apa = paste0(c("apa1_", "apa2_"), gene[1]), type = rep("apa", 2)))
       apas2 <- GRanges(seqnames = c(chr2[1], chr2[1]), strand = c(strand[1], strand[1]), 
                        ranges = IRanges(start = c(1100, 3305), width = c(1, 1)), 
                        DataFrame(gene = rep(gene2[1], 2), apa = paste0(c("apa1_", "apa2_"), gene2[1]), type = rep("apa", 2)))
       
       # We have two "identical" genes with the same two apa as the first test on chr1 and chr2
       features <- unlist(GRangesList(geneGr, geneGr2, apas1, apas2)) 
       # c() was complaining for disjoint chr (seqinfo)
       
       # FIXME if cigar=300M nothing was seen!
       rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(3301), 
                          cigar = "1M", strand = strand("+"))
       rd2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(3301), 
                          cigar = "1M", strand = strand("-"))
       rd3 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(3301), 
                          cigar = "1M", strand = strand("+"))
       rd4 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(3306), 
                          cigar = "1M", strand = strand("+"))
       rd5 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3306), 
                          cigar = "1M", strand = strand("+"))
       rd6 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3306), 
                          cigar = "1M", strand = strand("-"))
       rd7 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3306), 
                          cigar = "1M", strand = strand("+"))
       rd8 <- GAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3301), 
                          cigar = "1M", strand = strand("+"))
       
       reads <- suppressWarnings(list(c(rd1, rd2, rd3, rd4, rd5, rd6, rd7, rd8))) 
       # suppressWarnings for the same issue about different chr as before.
       # We have 3 reads on the PRE of the second choice APA and 1 on the post for gene A and the opposite for gene B.
       rds <- RoarDatasetMultipleAPA(reads, reads, features)
       rds <- countPrePost(rds)
       counts_geneA = matrix(c(0, 4, 0, 4, 3, 1, 3, 1), byrow = TRUE, ncol = 4)
       colnames(counts_geneA) = c("treatment_pre","treatment_post","control_pre","control_post")
       expect_equal(counts_geneA, assay(rds@roars[[1]],1))
       counts_geneB = matrix(c(0, 4, 0, 4, 1, 3, 1, 3), byrow = TRUE, ncol = 4)
       colnames(counts_geneB) = c("treatment_pre","treatment_post","control_pre","control_post")
       expect_equal(counts_geneB, assay(rds@roars[[2]],1))
       #rds <- computeRoars(rds)
       #rds <- computePvals(rds)
       #res <- countResults(rds)
    })

