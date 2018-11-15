#' BSAmpwrapper
#' Bunch of scripts that runs through the whole process
#' Might be prone to errors...


#' batchBSampseq
#'
#' Runs through the entire BSampseq script
#'
#' @param inputDir
#' @param outputDir
#'
batchBSampseq <- function(inputDir = ".", outputDir = ".", ...) {
  message("Finding manifest files...")
  man <- findManifest()
  if(length(man) == 0) return()
  man <- readManifest(man)
  message("Found and loaded manifest file!")
  txt <- findReports(inputDir)
  message("Found ", length(txt), " report files...")
  if(length(txt) == 0) return()
  if(!dir.exists(outputDir)) dir.create(outputDir)
  message("Processing reports...")
  batchprocessBSAmpSeq(txt, man, outputDir = outputDir, ...)
  fn <- list.files(outputDir)
  if(outputDir != ".") fn <- file.path(outputDir, fn)
  message("Creating standard curves...")
  gs <- guessStdBetas(fn)
  foo <- batchstdcurveBSampseq(gs)
}

#' batchstdcurveBSampseq
#' For creating multiple standard curves at once
#' Takes output from guessStdBetas() and will create a sample curve
#' for each ID
#' @param gs    output from guessStdBetas
#' @return nothing
batchstdcurveBSampseq <- function(gs) {
  gs <- gs[gs$IsStd,]
  stdlist <- unique(gs$ID)
  for (id in stdlist) {
    message("Making standard curve for sample ", id, "...")
    gssub <- gs[gs$ID %in% id,]
    foo <- stdcurveBSampseq(gssub$file, gssub$beta, id)

    xlname <- paste0(id, "_StdCurve.xlsx")
    std <- lapply(getSheetNames(xlname), read.xlsx, xlsxFile=xlname)
    names(std) <- getSheetNames(xlname)
    plotAllStdCurve(std, id)
  }
}
