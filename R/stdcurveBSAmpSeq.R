require(openxlsx)
require(stringr)
require(ggplot2)
require(gridExtra)
require(reshape2)

#' Get the names of amplicons from excel sheet names
.findSheetsWAmps <- function(file) {
  shts <- getSheetNames(file)
  shts <- shts[!sapply(shts, grepl, file)]
  if(length(shts) > 0) return(shts)
  # If this fails, just remove first sheet
  shts <- getSheetNames(file)
  shts <- shts[-1]
  return(shts)
}

#' Guess standards beta based on file name
#' P.S. This function does well but is not as smart as you are 
#' It will fail if the filenames are too obscure (or for other reasons)
#' PLEASE double check the results...
#' And when in doubt, just pick out the standards manually
#' @import readxl
#' @import openxlsx
#' @import stringr
#' @param filenames    names of files
#'                     can also take bsseq objects
#' @return data.frame with following columns
#'           ID: possible sample ID
#'           beta: beta value
#'           IsStd: is the sample a standard
#'           file: file path
guessStdBetas <- function(filenames) {
  if(class(filenames) == "BSseq") {
    myFilelist <- sampleNames(filenames)
  } else {
    myFilelist <- basename(filenames)
  }
  fnParse <- str_split_fixed(myFilelist, "[-_\\.]", 10)
  # Check if part of the file name has 0 and 100
  # samples <- unique(fnParse[,1])
  # fnIsStd <- sapply(samples, function(y) apply(fnParse[fnParse[,1] %in% y,], 2, function(x) any(x == 0) & any(x == 100)))
  fnIsStd <- apply(fnParse, 2, function(x) any(x == 0) & any(x == 100))
  fnIsStd[is.na(fnIsStd)] <- FALSE
  # SpIsStd <- apply(fnIsStd, 2, any)
  if(!any(fnIsStd) | sum(fnIsStd) > 1) {
    message("Unable to determine standards methylation value. ",
            "Methylation values are determined by a value embedded in ",
            "_ or .  For example, Standard_0_m.xlsx, ",
            "Standard_50_m.xlsx, Standard_100_m.xlsx")
    return(filenames)
  }
  # Get standard betas
  StdBeta <- fnParse[,fnIsStd]
  SampleIsStd <- grepl("^[0-9]+$", StdBeta)
  SampleIDs <- apply(fnParse[,!fnIsStd], 1, paste, collapse = "_")
  SampleIDs <- gsub("_targeted.*", "", SampleIDs)
  SampleIDs <- gsub("_*$", "", SampleIDs)
  # Find which samples have standards
  StdIDs <- sapply(unique(SampleIDs), function(x) all(SampleIsStd[SampleIDs %in% x]))
  StdIDs <- names(StdIDs)[StdIDs]
  # Save as bsseq or data.frame depending on input
  if(class(filenames) == "BSseq") {
    res <- filenames
    res$ID <- SampleIDs
    res$beta <- as.numeric(StdBeta)
  } else {
    res <- data.frame(ID = SampleIDs, beta = as.numeric(StdBeta),
                      IsStd = SampleIsStd, file = filenames,
                      stringsAsFactors=FALSE)
  }
  return(res)
}

#' stdcurveBSampseq
#' Generate standard curve for BS amplicon sequencing
#'
#' Outliers are defined as points which are +/- StdFitPercent
#' from the median beta value for each set of percent methylation
#' (i.e. median = 0.4, StdFitPercent = 0.4, acceptable range = 0-0.8)
#'
#' @import readxl
#' @import openxlsx
#'
#' @param bsseq         bsseq object
#'                      MUST have the following two pdata columns:
#'                       - ID  (sample ID, identical sample IDs will be plotted together)
#'                       - beta  (expected beta value)
#' @param regions       granges object of amplicons
#'                      MUST have the following column:
#'                       - genes  (gene associated with amplicon)
#' @param StdFitPercent (0.4) beta deviation from median to be classified as outlier
#' @param outputDir     (NULL) output folder for saving excels
#' @param saveExcel     (TRUE)  save excel file
#' @param return        (FALSE) returns value
#'
#' @return data.frame or excel file with standard curves
#'
#'
stdcurveBSampseq <- function(bsseq, regions, minCov = 5, StdFitPercent = 0.4,
                             outputDir = NULL, saveExcel = TRUE, return = FALSE) {
  # Order bsseq by increasing beta
  mySamples <- unique(bsseq$ID)
  mySamples <- mySamples[!is.na(mySamples)]
  
  # Do each standard individually
  for (mySample in mySamples) {
    # Get bsseq for subsample (remove 0 coverage regions)
    bs <- bsseq[,bsseq$ID %in% mySample]
    bs <- bs[rowSums(getCoverage(bs, type = "Cov")) != 0]
    # Filter by region coverage
    regCovs <- rowMeans(getCoverage(bs, regions = regions, type = "Cov", what = "perRegionAverage"), na.rm=T)
    reg <- regions[sapply(regCovs > minCov, isTRUE)]
    # Get methylation values
    M <- getMeth(bs, regions = reg, type = "raw")
    C <- getCoverage(bs, regions = reg, type = "Cov", what = "perBase")
    names(M) <- reg$genes
    names(C) <- reg$genes
    # Do this for each individual amplicon
    xlOut <- list()
    for (i in 1:length(M)) {
      # Reshape methylation
      gr2 <- granges(bs[bs %over% reg[i]])
      values(gr2) <- M[[i]]
      colnames(values(gr2)) <- sampleNames(bs)
      df <- melt(as.data.frame(gr2), 
                 id.vars = c("seqnames", "start", "end", "width", "strand"))
      names(df)[6:7] <- c("sample", "obs_beta")
      Stdbeta <- bs$beta
      names(Stdbeta) <- sampleNames(bs)
      df$exp_beta <- Stdbeta[df$sample]
      # Reshape coverage and merge with
      gr3 <-granges(bs[bs %over% reg[i]])
      values(gr3) <- C[[i]]
      df2 <- melt(as.data.frame(gr3), 
                 id.vars = c("seqnames", "start", "end", "width", "strand"))
      df$coverage <- df2$value
      
      # Identify outliers
      df$outlier <- FALSE
      CpGs <- unique(df$start)
      colMeds <- colMedians(M[[i]], na.rm = TRUE)
      for (iCpG in CpGs) {
        CpGbetas <- df$obs_beta[df$start == iCpG]
        if(any(abs(CpGbetas - colMeds) > StdFitPercent, na.rm = TRUE)) {
          df$outlier[df$start == iCpG] <- TRUE
        }
      }
      # Sort
      df <- df[order(df$exp_beta),]
      # Calculate standard curve
      PrimerRsq <- summary(lm(obs_beta ~ exp_beta, df[!df$outlier,], na.action=na.omit))$r.squared
      PrimerLM <- lm(obs_beta ~ exp_beta, df[!df$outlier,], na.action=na.omit)$coefficients
      df$variable <- c("R-squared", "Gradient", "Intercept", rep(NA, nrow(df)-3))
      df$value <- c(PrimerRsq, PrimerLM[2], PrimerLM[1], rep(NA, nrow(df)-3))
      # Save into list
      xlOut[[i]] <- df
    }
    names(xlOut) <- reg$genes
    
    # Write excel file
    outputfile <- paste0(mySample, "_StdCurve.xlsx")
    if(!is.null(outputDir)) outputfile <- file.path(outputDir, outputfile)
    if(saveExcel) write.xlsx(xlOut, file = outputfile, col.names = TRUE)
    if(return) return(xlOut)
  }

}




#' plotStdCurve
#' Plot standard curve for one amplicon
#'
#' @import ggplot2
#'
#' @param df          data.frame with data
#' @param line.yint   stdcurve y-intercept
#' @param line.slope  stdcurve slope
#' @param line.rsq    stdcurve r-squared
#' @param beta        name of observed beta column
#' @param expected    name of expected beta column
#' @param outlier     name of outlier column
#' @param coverage    logical: include coverage?
#' @param legend      logical: include legend?
#' @param title       title of plot
#' @param plotgraph   logical: plot the graph in R?
#' @param savefile    location of file to save or FALSE
#'
plotStdCurve <- function(df, line.yint, line.slope, line.rsq,
                         beta="obs_beta", expected="exp_beta", outlier="outlier",
                         coverage = TRUE,legend = TRUE, title,
                         plotgraph = FALSE, savefile = FALSE) {
  # Look for abline and rsq values if not defined
  if(all(missing(line.yint), missing(line.slope), missing(line.rsq))) {
    if(all(c("variable", "value") %in% names(df))) {
      myVars <- df$value
      names(myVars) <- df$variable
      line.yint <- myVars["Intercept"]
      line.slope <- myVars["Gradient"]
      line.rsq <- myVars["R-squared"]
    }
  }
  # Create base plot
  myPlot <- ggplot(df) +
    geom_point(aes_string(x=expected, y=beta, color=outlier),
               size = 2, na.rm = TRUE) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
    ylim(0,1)
  # Add title
  if(!missing(title)) {
    myPlot <- myPlot +
      labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  # Add best fit line
  if(!any(missing(line.yint), missing(line.slope))) {
    myPlot <- myPlot + geom_abline(intercept = line.yint, slope=line.slope)
  }
  # Add R-squared value
  if(!missing(line.rsq)) {
    myPlot <- myPlot + annotate("text", x=0, y=1, hjust=0, vjust=1, size=3, parse=TRUE,
                                label = paste("R^2 ==", round(line.rsq,3)))
  }
  if(coverage) {
    myPlot <- myPlot + annotate("text", x=0, y=0.937, hjust=0, vjust=1, size=3,
                                label = paste0("mean coverage = ", round(mean(df$coverage),1)))
  }
  # Exclude legend
  if(!legend) myPlot <- myPlot + theme(legend.position="none")
  # Show, save or return plot
  if(plotgraph) plot(myPlot)
  if(savefile != FALSE) ggsave(savefile, plot=myPlot, device="png",
                               width = 150, height = 100, units = "mm", dpi = 300)
  return(myPlot)
}


#' plotAllStdCurve
#' Wrapper for plotStdCurve
#' Generates plot for each amplicon, and an additional plot with all amplicons
#'
#' @import ggplot2
#' @import gridExtra
#'
#' @param dflist       output from stdcurveBSampseq
#'                     could be either list, or excel file
#' @param sampleID     name of sample displayed on plots
#' @param outputDir    output directory
#'
plotAllStdCurve <- function(dflist, sampleID="", outputDir=NULL, ...) {
  # If dflist is not a list but a character, try to load it as excel
  if(is.character(dflist) & all(grepl("xlsx$", dflist))) {
    dflist <- .read.xlsx.allsheets(dflist)
  }
  # Create a list of plots
  AmpPlots <- list()
  for (iAmp in 1:length(dflist)) {
    AmpDat <- dflist[[iAmp]]
    outputfile <- paste(sampleID, names(dflist[iAmp]), "StdCurve.png", sep="_")
    if(!is.null(outputDir)) outputfile <- file.path(outputDir, outputfile)
    plotStdCurve(AmpDat,
                 title = paste(sampleID, names(dflist)[iAmp], "Standard Curve"),
                 plotgraph = FALSE, savefile = outputfile, ...)
    AmpPlots[[iAmp]] <- plotStdCurve(AmpDat, title=names(dflist[iAmp]),
                                     plotgraph=FALSE, legend=FALSE)
  }
  # Create plot with all amplicons in one image
  names(AmpPlots) <- names(dflist)
  # Calculate height + width of plots
  iAmp <- round(sqrt(length(dflist)))
  iAmp[2] <- ceiling(length(dflist) / iAmp)
  # Create image and save file
  AmpGrobs <- arrangeGrob(grobs = AmpPlots, top=sampleID)

  outputfile <- paste(sampleID, "ALL_StdCurve.png", sep="_")
  if(!is.null(outputDir)) outputfile <- file.path(outputDir, outputfile)
  ggsave(outputfile, AmpGrobs,
         device="png", width = 100*iAmp[1], height = 75*iAmp[2],
         units = "mm", dpi=300)
}

#' Read all excel sheets
.read.xlsx.allsheets <- function(fn) {
  shts <- getSheetNames(fn)
  xl <- lapply(shts, read.xlsx , xlsxFile=fn)
  names(xl) <- shts
  return(xl)
}
