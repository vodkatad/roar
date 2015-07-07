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
#gene <- c("A", "B", NA, NA, NA, NA)
#type <- c("gene","gene","apa",  "apa", "apa", "apa")
#apa <- c(NA, NA, "apa1_A", "apa2_A", "apa1_B", "apa2_B")
#features <- GRanges(
#seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2","chr2")),
#strand = strand(rep("+", length(gene))),
#ranges = IRanges(
#start=c(1000, 2000, 1250, 1499, 2001, 2899),
#width=c(500, 900, 1,1,1,1)),
#DataFrame(gene, apa, type)
#)
#rd1 <- GAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))

