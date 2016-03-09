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
       
       # FIXME if cigar=300M nothing was seen! Due to the reads converted to single base!
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
       counts_geneA <- matrix(c(0, 4, 0, 4, 3, 1, 3, 1), byrow = TRUE, ncol = 4)
       colnames(counts_geneA) <- c("treatment_pre","treatment_post","control_pre","control_post")
       expect_equal(counts_geneA, assay(rds@roars[[1]],1))
       counts_geneB <- matrix(c(0, 4, 0, 4, 1, 3, 1, 3), byrow = TRUE, ncol = 4)
       colnames(counts_geneB) <- c("treatment_pre","treatment_post","control_pre","control_post")
       expect_equal(counts_geneB, assay(rds@roars[[2]],1))
       rds <- computeRoars(rds)
       rds <- computePvals(rds)
       res <- countResults(rds)
       mM <- c(-1,1,-1,-0.77777777)
       counts <- c(0, 3, 0, 1)
       lengths <- c(101, 6, 101, 6)
       countsResults <- data.frame(mM_treatment = mM, mM_control= mM, roar=rep(1,4), pval=rep(1,4),
                                   counts_treatment=counts, counts_control=counts, length=lengths)
       rownames(countsResults) <- c("A_apa1","A_apa2","B_apa1", "B_apa2")
       expect_equal(countsResults, res)
       fpkm <- c(0, 125000000, 0, 41666667)
       countsResults$length <- lengths
       countsResults$treatmentValue <- fpkm
       countsResults$controlValue <- fpkm
       countsResults$counts_control <- NULL
       countsResults$counts_treatment <- NULL
       res <- fpkmResults(rds)
       expect_equal(countsResults, res)
})

test_that("test_badIntegrationTestMinusStrand", 
   {
      chr <- rep("chr1", 4)
      strand <- rep("-", 4)
      gene <- rep("1", 4)
      type <- rep("gene", 4)
      apa <- rep(NA, 4)
      geneGr <- GRanges(seqnames = chr, strand = strand, ranges = IRanges(start = c(2, 4, 12, 19), width = c(1,4, 4, 3)), DataFrame(gene, apa, type))
      apas1 <- GRanges(seqnames = c(chr[1], chr[1], chr[1]), strand = c(strand[1], strand[1], strand[1]), 
                       ranges = IRanges(start = c(5, 7, 17), width = c(1, 1, 1)), 
                       DataFrame(gene = rep(gene[1], 3), apa = paste0(c("apa1_", "apa2_","apa3_"), gene[1]), type = rep("apa", 3)))
      features <- unlist(GRangesList(geneGr, apas1))
      
      r18_2 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(18), 
                         cigar = "1M", strand = strand("-"))
      r13_3 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(13), 
                         cigar = "4M", strand = strand("-"))
      r7_10 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(7), 
                         cigar = "1M", strand = strand("-"))
      r2_100 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2), 
                         cigar = "1M", strand = strand("-"))
      treatment_reads <- list(c(rep(r18_2,2), rep(r13_3,3), rep(r7_10, 10), rep(r2_100,100)))
      control_reads <- list(c(rep(r18_2,100), rep(r13_3,10), rep(r7_10, 3), rep(r2_100,2)))
      rds <- RoarDatasetMultipleAPA(treatment_reads, control_reads, features)
      rds <- countPrePost(rds)
      
      wanted_fragments_A <- GRanges(seqnames = rep(chr[1], 6), strand = rep(strand[1], 6), 
                                  ranges = IRanges(start = c(17,8,7,5,3,2), width = c(5, 9, 1, 2, 2, 1)), 
                                  DataFrame(length=c(3, 4, 1, 2, 1, 1)))
      expect_equal(wanted_fragments_A, rds@fragments[['1']])
      
      counts_geneA <- matrix(c(2, 113, 100, 15, 10, 100, 3, 2, 10, 100, 3, 2), byrow = TRUE, ncol = 4)
      colnames(counts_geneA) <- c("treatment_pre","treatment_post","control_pre","control_post")
      expect_equal(counts_geneA, assay(rds@roars[[1]],1))
      
      rds <- computeRoars(rds)
      rds <- computePvals(rds)
      res <- fpkmResults(rds)
      
      mM_t <- c(-0.946440938822624,-0.5921739130434782,-0.9307246376811594)
      mM_c <- c(19.579710144927535, 5.391304347826088, 0.1304347826086958)
      roars <- c(-0.04833784217524875, -0.10983870967741934, -7.135555555555548)
      fishers <- c(1.505474e-45, 0.01008237, 0.01008237)
      fpkm_c <- c(314465408.8050314, 28301886.79245283, 9433962.264150944)
      fpkm_t <- c(30303030.303030305, 454545454.54545456, 151515151.5151515)
      wanted_fpkmResults <- data.frame(mM_treatment = mM_t, mM_control = mM_c, roar = roars, pval = fishers, length=c(3,1,3),
                                       treatmentValue = fpkm_t, controlValue = fpkm_c)
      rownames(wanted_fpkmResults) <- c("1_apa3", "1_apa2", "1_apa1")
      expect_equal(wanted_fpkmResults, res, tolerance = .000002)
      
      # TODO ADD checks on filters to be safe from row order issues
})
