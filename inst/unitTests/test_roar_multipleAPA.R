# To manually perform single tests: source("inst/unitTests/test_roar.R")
# To perform all tests: 
# library("GenomicAlignments")
# BiocGenerics:::testPackage("roar")

test_RoarDatasetMultipleAPA_error_mismatch_gene_apa <- function() {
   gene <- c("A", "B", NA)
   type <- c("gene","gene","apa")
   apa <- c(NA, NA, "apa1_A")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr2", "chr1")),
      strand = strand(rep("+", length(gene))),
      ranges = IRanges(
         start=c(1000, 2000, 1300),
         width=c(500, 900, 1)),
      DataFrame(gene, apa, type)
   )
   
   rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   obs <- tryCatch(RoarDatasetMultipleAPA(list(c(rd1,rd1)), list(c(rd1,rd1)), features), 
                   error=function(e) e)
   checkEquals(obs$message, "All the genes in the gtf should have at least one apa")
}

#getApaGenesFractionsPlusStrand
test_getApaGenesFractionsPlusStrandBasic <- function() {
   chr <- rep("chr1", 3)
   strand <- rep("+", 3)
   gene <- rep("A", 3)
   type <- rep("gene", 3)
   apa <- rep(NA, 3)
   geneGr <- GRanges(
      seqnames = chr,
      strand = strand,
      ranges = IRanges(
         start=c(1000, 2000, 3300),
         width=c(500, 900, 10)),
      DataFrame(gene, apa, type)
   )
   apas <- GRanges(
      seqnames = c(chr[1], chr[1]),
      strand = c(strand[1], strand[1]),
      ranges = IRanges(
         start=c(1100, 3305),
         width=c(1,1)),
      DataFrame(gene=rep(gene[1],2), apa=paste0(c("apa1_", "apa2_"), gene[1]), type=rep("apa",2))
   )
   res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr[1], strand[1], gene[1])
   
   # Basic checks
   checkEqualsNumeric(length(res),2)
   checkEqualsNumeric(length(res[[1]]),4) # 4 fragments
   checkEqualsNumeric(length(res[[2]]),2) # 2 apaPrePost defs
   checkEqualsNumeric(res[[2]][[1]]@PREstart,1) # apaPrePost begin and end fragments
   checkEqualsNumeric(res[[2]][[1]]@PREend,1) # apaPrePost begin and end fragments
   checkEqualsNumeric(res[[2]][[2]]@PREstart,3) # apaPrePost begin and end fragments
   checkEqualsNumeric(res[[2]][[2]]@PREstart,3) # apaPrePost begin and end fragments
   checkEquals(res[[2]][[2]]@name,"apa2") # apaPrePost begin and end fragments
   checkEquals(res[[2]][[1]]@name,"apa1") # apaPrePost begin and end fragments
   # Fragments begin/ends and length (projected on exons):
   wanted <- GRanges(
      seqnames = rep(chr[1],4),
      strand = rep(strand[1],4),
      ranges = IRanges(
         start=c(1000, 1101,3300,3306),
         width=c(101,2199,6,4)),
      DataFrame(length=c(101,1299,6,4))
   )
   checkIdentical(wanted, res[[1]])
}

test_getApaGenesFractionsPlusStrandSingleIntron <- function() {
   chr <- "chr1"
   strand <- "+"
   gene <- "A"
   type <- "gene"
   apa <- NA
   geneGr <- GRanges(
      seqnames = rep(chr, 3),
      strand = rep(strand, 3),
      ranges = IRanges(
         start=c(1000, 2000, 3300),
         width=c(500, 900, 10)),
      DataFrame(rep(gene, 3), rep(apa, 3), rep(type, 3))
   )
   apas <- GRanges(
      seqnames = chr,
      strand = strand,
      ranges = IRanges(
         start=1510,
         width=1),
      DataFrame(gene=gene, apa=paste0("apa1_", gene), type="apa")
   )
   res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr, strand, gene)
   
   # Basic checks
   checkEqualsNumeric(length(res[[1]]),3) # 4 fragments
   checkEqualsNumeric(length(res[[2]]),1) # 2 apaPrePost defs
   checkEqualsNumeric(res[[2]][[1]]@PREstart,1) # apaPrePost begin and end fragments
   checkEqualsNumeric(res[[2]][[1]]@PREend,1) # apaPrePost begin and end fragments
   checkEquals(res[[2]][[1]]@name,"apa1") # apaPrePost begin and end fragments
   
   # Fragments begin/ends and length (projected on exons):
   wanted <- GRanges(
      seqnames = rep(chr,3),
      strand = rep(strand,3),
      ranges = IRanges(
         start=c(1000,1511,3300),
         width=c(511,1789,10)),
      DataFrame(length=c(500,900,10))
   )
   checkIdentical(wanted, res[[1]])
   
   ## I would have expected only 2 fragments but the insertion of a
   ## fake APA representing the gene end (to exploit the already
   ## written logic) will always put a fragment == to the last exon.
   ## XXX TODO understand if this causes problems or slows down the
   ## execution too much (I don't think so).
   
   # Single APA in exon (and in intron for that matter) should be extensively tested
   # in the comparison with the fake multiple gtf with a single APA (most "distal")
   # for each gene.
}

test_getApaGenesFractionsPlusStrandIntronExon <- function() {
   chr <- "chr1"
   strand <- "+"
   gene <- "A"
   type <- "gene"
   apa <- NA
   geneGr <- GRanges(
      seqnames = rep(chr, 3),
      strand = rep(strand, 3),
      ranges = IRanges(
         start=c(1000, 2000, 3300),
         width=c(500, 900, 10)),
      DataFrame(rep(gene, 3), rep(apa, 3), rep(type, 3))
   )
   apas <- GRanges(
      seqnames = rep(chr,2),
      strand = rep(strand,2),
      ranges = IRanges(
         start=c(2010,2901),
         width=c(1, 1)),
      DataFrame(gene=rep(gene,2), apa=paste0(c("apa1_","apa2_"), gene), type=rep("apa",2))
   )
   res <- getApaGenesFractionsPlusStrand(geneGr, apas, chr, strand, gene)
   
   # Basic checks
   checkEqualsNumeric(length(res[[1]]),4)
   checkEqualsNumeric(length(res[[2]]),2)
   checkEqualsNumeric(res[[2]][[1]]@PREstart,1)
   checkEqualsNumeric(res[[2]][[1]]@PREend,1)
   checkEqualsNumeric(res[[2]][[2]]@PREstart,1)
   checkEqualsNumeric(res[[2]][[2]]@PREend,2)

   
   # Fragments begin/ends and length (projected on exons):
   wanted <- GRanges(
      seqnames = rep(chr[1],4),
      strand = rep(strand[1],4),
      ranges = IRanges(
         start=c(2000,2011,2902,3300),
         width=c(11,891,398,10)),
      DataFrame(length=c(11,889,0,10))
   )
   checkIdentical(wanted, res[[1]])
}

# Situations to test:
# single apa? In intron, in exon. Multiple apas in same exon, in same intron. Last/first?
# minus strand
# counts 
# A simple sample gene with 3 exons is sufficient? Are there extreme cases for intron/exon structure? 
# Apa on the boundaries.

# unordered calls gets fixed or not
# paired pvals
# multiple samples  generateRoarsMultipleBam
# XXX decide what to do about lengths == 0, remove? Should be filtered a priori (better/simpler choice maybe).