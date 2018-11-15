

#' 
#' 
#' @param dir
loadbsseqBSAmpSeq <- function(dir=".", gr=NULL) {
  fn <- list.files(dir, "cov.gz")
  fn <- file.path(dir, fn)
  sn <- gsub("_S.*", "", basename(fn))
  bsseq <- read.bismark(fn, sn, strandCollapse = FALSE)
  bsseq <- addStrandFromBSGRanges(bsseq, gr)
  return(bsseq)
}

#' Build granges object with all CpG loci
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @return granges object with all CpG loci
buildEmptyBSseqGranges <- function (genome="hg19") {
  require(GenomicRanges)
  require(paste0("BSgenome.Hsapiens.UCSC.", genome), character.only = TRUE)
  chrs <- names(Hsapiens)[1:24]
  cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
  cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))
  seqinfo(cpgr) <- seqinfo(Hsapiens)[seqlevels(cpgr)]
  return(cpgr)
}

#' Adds strand information into a bsseq object using a granges object
#' containing all CpG loci (i.e. granges(BSseq-object))
#' If no granges object is provided, it will use
#' "buildEmptyBSseqGranges" function to generate one
#' @param bs      bsseq object
#' @param gr      granges object
#' @param genome  genome (if no granges is provided)
#' @return bsseq object with strand information
addStrandFromBSGRanges <- function(bs, gr=NULL, genome="hg19") {
  if(is.null(gr)) {
    message("Generating strand information from reference genome may take some time...")
    gr <- buildEmptyBSseqGranges(genome)
  }
  iscpg <- from(findOverlaps(granges(bs), gr))
  bs <- bs[iscpg,]
  # Find strand
  isplus <- findOverlaps(granges(bs), gr, type = "start")
  # isminus <- findOverlaps(granges(bs), gr, type = "end")
  strand(bs) <- Rle(ifelse(1:length(bs) %in% from(isplus), "+", "-"))
  return(bs)
}

#' bsseq2ampseqreport
#' Generates CX_report.txt file similar to basespace targeted
#' amplicon sequencing format from a bsseq object
#' @param bsseq     bsseq object
#' @param return    (FALSE) returns data.frame
#' @return 
bsseq2ampseqreport <- function(bsseq, return=FALSE) {
  sn <- sampleNames(bsseq)
  for (i in sn) {
    bsseqS <- bsseq[,i]
    # Remove loci with no coverage
    M <- getMeth(bsseqS, type="raw")
    bsseqS <- bsseqS[!as.vector(is.nan(M)),]
    # Calculate M and U values
    M <- getCoverage(bsseqS, type="M")
    U <- getCoverage(bsseqS, type="Cov") - M
    # out <- data.frame(data.frame(granges(bsseqS))[,c(1,2,5)], M, U, "CG", "CG?")
    out <- data.frame(data.frame(granges(bsseqS))[,c(1,2,5)], M, U, "CG")
    write.table(out, paste0(i, ".targeted.CX_report.txt"), sep="\t",
                quote = F, row.names = F, col.names = F)
    if(return) return(out)
  }
}