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
test_getApaGenesFractionsPlusStrand <- function() {
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
   
   # XXX fare xls preciso length e inizi/fini attesi e confronto per i fragments
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

