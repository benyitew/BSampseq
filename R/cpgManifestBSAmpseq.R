

#' Creates a CpG manifest from BSAmpSeq standard curve file
#' @param std      standard curve excel file
#' @param save     either FALSE or the file name to save excel as
#' @param return   logical - return data.frame
createCpgManifest <- function(std, save=FALSE, return=FALSE) {
  if(class(std) == "character") std <- .read.xlsx.allsheets(std)
  std <- mapply(function(x,y) data.frame(Name=y, x), x=std, y=names(std), SIMPLIFY = F)
  std <- Reduce(rbind, std)
  std$Coord <- paste0(substr(std$Chrm,4,5), "-", std$Position)
  cpgs <- std[!duplicated(std$Coord),c(1:4,10)]
  if(!isFALSE(save)) write.xlsx(cpgs, save)
  return(cpgs)
}

#' Converts CpG manifest file to granges object
#' @param man     manifest file
#' @return granges object
CpGManifest2GRanges <- function(man) {
  names(man)[2:4] <- c("chr", "start", "strand") 
  man$end <- man$start
  data.frame2GRanges(man, keepColumns = T)
}


#' getMeth or getCovg on bsseq object based on ranges given in CpG manifest
#' @param bs     bsseq file
#' @param man    manifest file
#' @param exclude_outliers  logical - exclude outliers
#' @param show_coords       logical - include coordinates in output
#' @return data.frame 
getMethFromCpgManifest <- function(bs, man, exclude_outliers = TRUE, show_coords = FALSE) {
  gr <- CpGManifest2GRanges(man)
  if(isTRUE(exclude_outliers)) gr <- gr[!gr$Outlier]
  M <- getMeth(bs, regions = gr, type="raw", what = "perRegion")
  M <- data.frame(name = gr$Name, M)
  if(isTRUE(show_coords)) M <- cbind(data.frame(gr)[,c(1,2,5)], M)
  return(M)
}
getCovgFromCpgManifest <- function(bs, man, exclude_outliers = TRUE, show_coords = FALSE) {
  gr <- CpGManifest2GRanges(man)
  if(isTRUE(exclude_outliers)) gr <- gr[!gr$Outlier]
  C <- getCoverage(bs, regions = gr, type="Cov", what = "perRegionTotal")
  C <- data.frame(name = gr$Name, C)
  if(isTRUE(show_coords)) C <- cbind(data.frame(gr)[,c(1,2,5)], C)
  return(C)
}

#' Summarizes output from getMeth/CovgFromCpGManifest by amplicon name
#' @param df    data.frame containing meth or covg values
summarizeAmpliconByName <- function(df) {
  v <- df[,-which(names(df) %in% c("seqnames", "start", "strand", "name"))]
  res <- sapply(unique(df$name), function(x) colMeans(v[df$name %in% x,], na.rm = T))
  data.frame(name=unique(df$name), t(res))
}
