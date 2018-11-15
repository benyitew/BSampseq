require(openxlsx)
require(stringr)
require(ggplot2)
require(gridExtra)

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
#' @return data.frame with following columns
#'           ID: possible sample ID
#'           beta: beta value
#'           IsStd: is the sample a standard
#'           file: file path
guessStdBetas <- function(filenames) {
  myFilelist <- basename(filenames)
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
    return(NULL)
  }
  # Get standard betas
  StdBeta <- fnParse[,fnIsStd]
  SampleIsStd <- grepl("^[0-9]+$", StdBeta)
  SampleIDs <- apply(fnParse[,!fnIsStd], 1, paste, collapse = "_")
  SampleIDs <- gsub("_targeted.*", "", SampleIDs)
  # Find which samples have standards
  StdIDs <- sapply(unique(SampleIDs), function(x) all(SampleIsStd[SampleIDs %in% x]))
  StdIDs <- names(StdIDs)[StdIDs]

  res <- data.frame(ID = SampleIDs, beta = as.numeric(StdBeta),
                    IsStd = SampleIsStd, file = filenames,
                    stringsAsFactors=FALSE)
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
#' @param filenames     processed excel files (use processBSampseq)
#' @param betas         beta values of files
#' @param name          name of sample
#' @param StdFitPercent (0.4) beta deviation from median to be classified as outlier
#' @param outputDir     (NULL) output folder for saving excels
#' @param saveExcel     (TRUE)  save excel file
#' @param return        (FALSE) returns value
#'
#' @return data.frame or excel file with standard curves
#'
#'
stdcurveBSampseq <- function(filenames, betas, name, StdFitPercent = 0.4,
                             outputDir = NULL, saveExcel = TRUE, return = FALSE) {
  # Order files by increasing beta
  myFileList <- filenames[order(betas)]
  StdBeta <- betas[order(betas)]

  # Find common amplicons in all standards file
  mySheets <- lapply(myFileList, .findSheetsWAmps)
  mySheets.common <- sapply(mySheets, function(x) mySheets[[1]] %in% x)
  if(class(mySheets.common) == "logical") {
    mySheets.common <- all(mySheets.common)
  } else {
    mySheets.common <- apply(mySheets.common, 1, all)
  }
  if(!all(mySheets.common)) message("Some amplicons were not found in all standards, only common amplicons were used.")
  mySheets <- mySheets[[1]][mySheets.common]

  # Read data for all standards
  Stds <- lapply(myFileList, function(x) lapply(mySheets, function(y) read.xlsx(x, sheet = y, colNames = TRUE)))
  # If no correct header, but has over 8 columns
  if(!"Beta" %in% names(Stds[[1]][[1]])) {
    Stds <- lapply(myFileList, function(x) lapply(mySheets, function(y) read.xlsx(x, sheet = y, colNames = FALSE)))
    if(all(sapply(Stds[[1]], ncol) >= 7)) {
      Stds <- lapply(Stds, function(x) lapply(x, function(y) setNames(y, c("Chrm","Position","Strand","Methyl Read","UnMethyl Reads","Context","Seq","Beta"))))
    } else {
      stop("Error reading standards, make sure the columns are properly labeled.")
    }
  }
  names(Stds) <- StdBeta
  Stds <- lapply(Stds, function(x) setNames(x, mySheets))

  # Identify only common columns and trim off excess columns
  myCols <- lapply(Stds, function(x) sapply(x, function(y) names(y)))
  myCols <- sapply(unique(as.vector(myCols[[1]])), function(z)
    all(sapply(myCols, function(x) all(apply(x, 2, function(y) z %in% y))))
  )
  myCols <- names(myCols[myCols])
  myCols <- myCols[myCols %in% c("Chrm","Position","Strand","Meth_Reads","UnMeth_Reads","Context","Seq","Beta")]
  Stds <- lapply(Stds, function(x) lapply(x, function(y) y[,myCols]))
  if(!"Beta" %in% myCols) stop("Beta column not found in one or more files.")

  # Add percent methylation value
  for (i in 1:length(Stds)) {
    Stds[[i]] <- lapply(Stds[[i]], function(x) data.frame(x, Exp_Beta=StdBeta[i]))
  }


  # Combine all standards for each amplicon
  AmpStds <- list()
  for (i in 1:length(mySheets)) {
    AmpStds[[mySheets[i]]] <- do.call("rbind", sapply(Stds, function(x) x[i]))
  }

  # Identify outliers for each amplicon
  for(iAmp in 1:length(AmpStds)) {
    AmpDat <- AmpStds[[iAmp]]
    # Create a list of CpGs, using position as name
    PrimerGood <- rep(TRUE, length(unique(AmpDat$Position)))
    names(PrimerGood) <- unique(AmpDat$Position)
    # Go through each CpG and check if the median is greater than StdFitPercent
    for (i in unique(AmpDat$Exp_Beta)) {
      CpGNums <- match(AmpDat$Position[AmpDat$Exp_Beta %in% i], names(PrimerGood))
      OneBeta <- AmpDat$Beta[AmpDat$Exp_Beta %in% i]
      PrimerGood[CpGNums] <- !abs(OneBeta - median(OneBeta)) > StdFitPercent & PrimerGood[CpGNums]
    }
    AmpDat$Outlier <- !PrimerGood[match(AmpDat$Position, names(PrimerGood))]
    AmpStds[[iAmp]] <- AmpDat
  }
  # Check if there are no values in amplicon
  AmpLMStds <- AmpStds
  AmpYesReads <- sapply(AmpStds, function(x) any(!x$Outlier))
  AmpYesReads[is.na(AmpYesReads)] <- FALSE
  if(any(!AmpYesReads)) {
    AmpLMStds <- AmpStds[AmpYesReads]
    print(paste("No reads found in the following amplicons, no standard curve were made:",
                paste(names(which(!AmpYesReads)), collapse=", ")))

  }
  # Calculate linear model and r-square (exclude outliers)
  PrimerRsq <- sapply(AmpLMStds, function(x) summary(lm(Beta ~ Exp_Beta, x[!x$Outlier,], na.action=na.omit))$r.squared)
  PrimerLM <- sapply(AmpLMStds, function(x) lm(Beta ~ Exp_Beta, x[!x$Outlier,], na.action=na.omit)$coefficients)
  # Add these numbers to matrices for export
  for (i in names(AmpLMStds)) {
    AmpStds[[i]]$variable <- c("R-squared", "Gradient", "Intercept", rep(NA, nrow(AmpStds[[i]])-3))
    AmpStds[[i]]$value <- c(PrimerRsq[i], PrimerLM[2,i], PrimerLM[1,i], rep(NA, nrow(AmpStds[[i]])-3))
  }

  # Write excel file
  outputfile <- paste(name, "_StdCurve.xlsx", sep="")
  # outputfile <- gsub(pattern = ".xlsx$", replacement = "_StdCurve.xlsx", myFileList[length(myFileList)])
  # outputfile <- gsub(pattern = "[\\_\\.]100", replacement = "", outputfile)
  if(!is.null(outputDir)) outputfile <- file.path(outputDir, outputfile)
  if(saveExcel) write.xlsx(AmpStds, file = outputfile, col.names = TRUE)
  if(return) return(AmpStds)
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
#' @param legend      logical: include legend?
#' @param title       title of plot
#' @param plotgraph   logical: plot the graph in R?
#' @param savefile    location of file to save or FALSE
#'
plotStdCurve <- function(df, line.yint, line.slope, line.rsq,
                         beta="Beta", expected="Exp_Beta", outlier="Outlier",
                         legend = TRUE, title,
                         plotgraph=FALSE, savefile=FALSE) {
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
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"))
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
plotAllStdCurve <- function(dflist, sampleID="", outputDir=NULL) {
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
                 plotgraph = FALSE, savefile = outputfile)
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
