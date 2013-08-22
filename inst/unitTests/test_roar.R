# To manually perform single tests: source("inst/unitTests/test_roar.R")

test_test <- function() {
   checkEquals(2, 2)   
   #checkTrue(is.na(divideBy(4, 0)))
   #checkEqualsNumeric(divideBy(4, 1.2345), 3.24, tolerance=1.0e-4)
}

test_countPrePost_singleSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", "D_POST", "E_PRE", "E_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2", "chr2", "chr2",
                       "chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000),
         width=c(500, 900, 500, 300, 600, 300, 500, 900)),
      DataFrame(gene_id)
   )
   rd1 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(c(rd1,rd2)), list(c(rd3,rd4)), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds,1)[1,1], 1)
   checkEqualsNumeric(assay(rds,1)[1,2], 1)
   checkEqualsNumeric(assay(rds,1)[2,3], 2)
}

test_countPrePost_mulSamples <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST", "D_PRE", "D_POST", "E_PRE", "E_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2", "chr2", "chr2",
                       "chr1", "chr1")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000),
         width=c(500, 900, 500, 300, 600, 300, 500, 900)),
      DataFrame(gene_id)
   )
   rd1 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(2000), cigar = "300M", strand = strand("+"))
   rd3 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(rd1,rd2), list(rd3,rd4), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds@countsRight,1)[1,1], 1)
   checkEqualsNumeric(assay(rds@countsRight,1)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsLeft,2)[2,1], 1)
}

test_countPrePost_preferPOST <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   # Had to change from Rle() to factor for seqnames otherwise there was a warning about the list for the left
   # alignments about having only chr1 for some reads and chr2 for the others.
   rd1 <- GappedAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
   rd2 <- GappedAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "3000M", strand = strand("+"))
   rd3 <- GappedAlignments("a", seqnames = factor("chr2", levels=c("chr1","chr2")), pos = as.integer(2800), cigar = "300M", strand = strand("+"))
   rd4 <- GappedAlignments("a", seqnames = factor("chr2", levels=c("chr1","chr2")), pos = as.integer(3500), cigar = "300M", strand = strand("+"))
   rd5 <- GappedAlignments("a", seqnames = factor("chr1", levels=c("chr1","chr2")), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(c(rd1,rd2)), list(c(rd3,rd4,rd5)), features)
   rds <- countPrePost(rds, FALSE)
   checkEqualsNumeric(assay(rds,1)[1,1], 0)
   checkEqualsNumeric(assay(rds,1)[1,2], 2)
   checkEqualsNumeric(assay(rds,1)[1,3], 1)
   checkEqualsNumeric(assay(rds,1)[1,4], 0)
   checkEqualsNumeric(assay(rds,1)[2,3], 1)
   checkEqualsNumeric(assay(rds,1)[2,4], 1)
   checkEqualsNumeric(assay(rds,1)[2,1], 0)
   checkEqualsNumeric(assay(rds,1)[2,2], 0)
}

test_countPrePost_stranded <- function() {
   gene_id <- c("A_PRE", "A_POST", "B_PRE", "B_POST")
   features <- GRanges(
      seqnames = Rle(c("chr1", "chr1", "chr2", "chr2")),
      strand = strand(rep("+", length(gene_id))),
      ranges = IRanges(
         start=c(1000, 2000, 3000, 3600),
         width=c(1000, 900, 600, 300)),
      DataFrame(gene_id)
   )
   rd1 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("+"))
   rd2 <- GappedAlignments("a", seqnames = Rle("chr1"), pos = as.integer(1000), cigar = "300M", strand = strand("-"))
   rd3 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3000), cigar = "300M", strand = strand("+"))
   rd4 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("-"))
   rd5 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   rd6 <- GappedAlignments("a", seqnames = Rle("chr2"), pos = as.integer(3200), cigar = "300M", strand = strand("+"))
   
   rds <- RoarDataset(list(rd1,rd2), list(rd3,rd4,c(rd5,rd6)), features)
   rds <- countPrePost(rds, TRUE)
   checkEqualsNumeric(assay(rds@countsRight,1)[1,1], 1)
   checkEqualsNumeric(assay(rds@countsRight,1)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsRight,2)[1,1], 0)
   checkEqualsNumeric(assay(rds@countsRight,2)[1,2], 0)
   checkEqualsNumeric(assay(rds@countsLeft,1)[2,1], 1)
   checkEqualsNumeric(assay(rds@countsLeft,2)[2,1], 0)
   checkEqualsNumeric(assay(rds@countsLeft,3)[2,1], 2)
}