require(openxlsx)

#' Description:
#'  - Identify amplicons of interest from targeted alignment coordinates file
#'   (The first xlsx file with the name "alignment" OR "manifest" will be used)
#'  - Extract amplicons from targeted report files
#'  - Save results as excel file
#'
#' Inputs ## Targeted report text files with the following columns:
#' Chrm, Position, Strand, Methyl Read, UnMethyl Reads, Context, Seq
#' There should be no column headers!
#'
#' Targeted alignment coordinates / manifest file must have the following columns:
#' Name, Strand, Chromosome, Amplicon Start, Amplicon End
#'
#' Outputs ## Combined excel file for each amplicon
#'

#' findReports
#' @param dir        directory to look for reports
#' @param fileregex  regex string to look for reports
#' @return  vector containing (relative) report filepath
findReports <- function(dir = ".", fileregex = "txt$") {
  myFileList <- list.files(dir, pattern = fileregex)
  if(dir != ".") myFileList <- file.path(dir, myFileList)
  return(myFileList)
}

#' findManifest
#' @param dir        directory to look for manifest
#' @param fileregex  regex string to look for manifest
#' @return  string containing (relative) manifest filepath
findManifest <- function(dir = ".", fileregex = ".*([Mm]anifest|[Aa]lignment).*xlsx$") {
  # Find manifest file
  myManifest <- list.files(dir, fileregex)
  if(length(myManifest) > 1) {
    message("Multiple manifest files found, the following manifest file was used: ", myManifest[1])
    myManifest <- myManifest[1]
  }
  # Compile file path if necessary
  if (length(myManifest) == 1) {
    if (dir != ".") myManifest <- file.path(dir, myManifest)
  }
  if (length(myManifest) < 1) message("No manifest files found. Please specify manifest file location.")
  return(myManifest)
}

#' readManifest
#' Reads a manifest file
#' Manifest files MUST contain the following columns:
#'   - Name, Strand, Chromosome, Amplicon Start, Amplicon End
#' Also does the following conversions:
#'   - Strand: "plus", "minus" to "+", "-"
#'   - Chromosome: "1", "21" to "chr1", "chr21"
#' @import openxlsx
#' @param manifest   Can be either:
#'                    string - path to manifest file
#'                    data.frame - actual manifest file
#' @return formatted manifest data.frame
readManifest <- function (myManifest) {
  # Read excel if file path is given
  if(is.character(myManifest)) myManifest <- read.xlsx(myManifest)
  # Find the required columns and arrange
  RequiredColumns <- c("Name", "Strand", "Chromosome", "Start", "End")
  whichCols <- sapply(RequiredColumns, function(x) which(grepl(x, names(myManifest), ignore.case = TRUE)))
  if(is.list(whichCols)) {
    missingCol <- names(which(sapply(whichCols, length) == 0))
    message("Unable to find the following columns: ", paste(missingCol, collapse = ", "))
    return()
  }
  myManifest <- myManifest[,whichCols]
  names(myManifest) <- RequiredColumns

  # Modify columns as needed
  myManifest$Strand[myManifest$Strand %in% "plus"] <- "+"
  myManifest$Strand[myManifest$Strand %in% "minus"] <- "-"
  if(is.numeric(myManifest$Chromosome)) {
    myManifest$Chromosome <- paste("chr", myManifest$Chromosome, sep="")
  }
  return(myManifest)
}

#' processBSampseq
#' Processes one bisulfite-amplicon sequencing file
#'
#' @import openxlsx
#' @param filename      file path
#' @param manifest      data.frame containing manifest
#' @param minReads      (10) minimum number of reads required
#' @param outputDir     output folder for saving excels
#' @param saveExcel     (TRUE)  save excel file
#' @param return        (FALSE) returns value
#' @return excel and/or value depending on arguments
processBSampseq <- function(filename, manifest, minReads = 10,
                            outputDir = NULL, saveExcel = TRUE, return = FALSE) {
  # Read file and check if file looks valid
  methylreads <- read.delim(filename, header = FALSE)
  if(!"integer" %in% sapply(methylreads, class)) methylreads <- read.delim(filename, header = TRUE)
  if(!"integer" %in% sapply(methylreads, class)) {
    message("Unable to process file:", filename)
    return(NULL)
  }
  # Get sample ID from file name & generate output file
  sampleID <- gsub("\\_[A-Za-z0-9]+\\..*", "", x = basename(filename))
  # Format input
  names(methylreads) <- c("Chrm","Position","Strand","Meth_Reads","UnMeth_Reads","Context","Seq")
  methylreads <- methylreads[,1:7]
  # Create list for tabbed excel file, starting with all reads
  myOutput <- list()
  myOutput[[sampleID]]<- methylreads
  # Keep only reads where (context == CG)
  methylreads <- methylreads[methylreads$Context %in% "CG",]
  # Find all regions in manifest file
  for (iAmp in 1:nrow(manifest)) {
    # Filter by: chromosome, strand and position
    mHits <- methylreads[methylreads$Chrm == manifest$Chromosome[iAmp],]
    mHits <- mHits[mHits$Strand == manifest$Strand[iAmp],]
    mHits <- mHits[mHits$Position >= manifest$Start[iAmp],]
    mHits <- mHits[mHits$Position <= manifest$End[iAmp],]
    # Calculate Beta
    mHits$Beta <- mHits$Meth_Reads/(mHits$Meth_Reads + mHits$UnMeth_Reads)
    # If there are hits after filtering, save sheet
    if (nrow(mHits) > 0 & any(mean(mHits$Meth_Reads) >= minReads,
                              mean(mHits$UnMeth_Reads) >= minReads)) {
      myOutput[[as.character(manifest$Name[iAmp])]] <- mHits
    }
  }
  # Save as excel or return nvalue
  if(saveExcel) {
    outputfile <- gsub(pattern = "txt$", replacement = "xlsx", basename(filename))
    if(!is.null(outputDir)) outputfile <- file.path(outputDir, outputfile)
    # write.xlsx(myOutput, file = outputfile, col.names = TRUE)
    write.xlsx(myOutput, file = outputfile)
  }
  if(return) return(myOutput)
}
#' batchprocessBSAmpSeq
#' Loops processBSampseq() over multiple files
#' @param files    vector of files to be processed
#' @param ...      arguments passed to processBSampseq
batchprocessBSAmpSeq <- function(files, ...) {
  foo <- lapply(files, processBSampseq, ...)
}

