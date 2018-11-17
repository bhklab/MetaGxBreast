#' Function to load breast cancer expression sets from the Experiment Hub
#'
#' This function returns breast cancer datasets from the hub and a vector of patients from the datasets that are most likely duplicates
#' @param loadString a character vector specifying which data will be loaded. The default is "majority", which loads in 37 of the 39 datasets.
#' The other option is to provide a character vecotr of the names of the datasets to load. The metabric and tcga datasets areloaded separately as they are 
#' very large and doing so will help prevent memory allocation errors for R windows. Furthermore, these datasets are so large that they
#' dominate statistical analyses so it is best that they are analyzed separate of the 37 smaller datasets loaded with the string majority
#' @param removeDuplicates remove patients with a Spearman correlation greater than or equal to 0.98 with other patient expression profiles (default TRUE)
#' @param quantileCutoff A nueric between 0 and 1 specifying to remove genes with standard deviation below the required quantile (default 0)
#' @param rescale apply centering and scaling to the expression sets (default FALSE)
#' @param minNumberGenes an integer specifying to remove expression sets with less genes than this number (default 0)
#' @param minNumberEvents an integer specifying how man survival events must be in the dataset to keep the dataset (default 0)
#' @param minSampleSize an integer specifying the minimum number of patients required in an eset (default 0)
#' @param removeRetracted remove datasets from retracted papers (default TRUE, currently just PMID17290060 dataset)
#' @param removeSubsets remove datasets that are a subset of other datasets (defeault TRUE, currently just PMID19318476)
#' @param keepCommonOnly remove probes not common to all datasets (default FALSE)
#' @param imputeMissing remove patients from datasets with missing expression values
#' @return a list with 2 elements. The First element named esets contains the datasets. The second element named duplicates contains
#' a vector with patient IDs for the duplicate patients (those with  Spearman correlation greater than or equal to 0.98 with other patient expression profiles).
#' @export
#' @importFrom Biobase esApply featureNames sampleNames exprs pData experimentData
#' @importFrom lattice levelplot
#' @importFrom impute impute.knn
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @importFrom stats complete.cases sd quantile
#' @examples
#'
#' #Use the default loadString = "majority" if you want the 37 smaller datasets
#' esetsAndDups = loadBreastEsets(loadString = c("CAL", "DFHCC", "DFHCC2", "DFHCC3", "DUKE", "DUKE2", "EMC2"))


loadBreastEsets = function(loadString = "majority", removeDuplicates = TRUE, quantileCutoff = 0, rescale = FALSE, minNumberGenes = 0,
                           minNumberEvents = 0, minSampleSize = 0, removeRetracted = TRUE, removeSubsets = TRUE,
                           keepCommonOnly = FALSE, imputeMissing = FALSE)
{
  duplicates = NULL
  #if(getRversion() >= "2.15.1")  utils::globalVariables(c("."), add = F)
  ## -----------------------------------------------------------------------------
  ## needed functions
  ## -----------------------------------------------------------------------------
  filterQuantile <- function(object, q){
    if (!identical(q >=0 && q < 1, TRUE))
      stop("require 0 <= q < 1")
    if (!identical(class(object) == "ExpressionSet", TRUE))
      stop("object must be an ExpressionSet")
    geneSd <- Biobase::esApply(object,1,sd, na.rm=TRUE)
    gene.quantile <- stats::quantile(geneSd, probs=q)
    actual.makescutoff <- sum(geneSd < gene.quantile) / length(geneSd)
    ##make sure the correct number of genes are getting filtered:
    if (abs(q - actual.makescutoff) > 0.01){
      stop("Not scaling this object, likely pre-scaled.")
    }else{
      object <- object[geneSd > gene.quantile, ]
    }
    return(object)
  }
  ##recursive intersect function
  intersectMany <- function(lst){
    ## Find the intersection of multiple vectors stored as elements of a
    ## list, through a tail-recursive function.
    if (length(lst)==2){
      return(intersect(lst[[1]],lst[[2]]))
    }else{
      return(intersectMany(c(list(intersect(lst[[1]],lst[[2]])),lst[seq(-1, -2)])))
    }
  }
  
  ##Split out non-specific probe sets
  expandProbesets <- function (eset, sep = "///"){
    x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
    eset <- eset[order(vapply(x, length, numeric(1))), ]
    x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
    idx <- unlist(vapply(x, function(i) rep(i, length(x)), character(length(x))))
    xx <- !duplicated(unlist(x))
    idx <- idx[xx]
    x <- unlist(x)[xx]
    eset <- eset[idx, ]
    Biobase::featureNames(eset) <- x
    eset
  }
  
  ## -----------------------------------------------------------------------------
  ##load the esets
  ## -----------------------------------------------------------------------------
  
  hub = ExperimentHub::ExperimentHub()
  #AnnotationHub::possibleDates(hub)
  #ovarianData = query(hub, c("MetaGxOvarian", "ExpressionSet"))
  breastData = query(hub, c("MetaGxBreast", "ExpressionSet"))
  tcgaInd = which(grepl("TCGA", breastData$title))
  metabricInd = which(grepl("METABRIC", breastData$title))
  if(length(loadString) == 1)
  {
    if(loadString == "majority"){
      breastData = breastData[seq_len(length(breastData))[-c(metabricInd, tcgaInd)]]
    }else{
      stop("loadString needs to be one of majority, metabric, tcga, or a character vector of datasets to load from the hub")
    }
  }else{
    keepIndVec = c()
    for(i in seq_len(length(loadString)))
    {
      keepInd = which(grepl(paste0(loadString[i], "_"), paste0(breastData$title, "_")))
      if(length(keepInd) == 0){
        stop(paste(loadString[i], "could not be found in the MetaGxBreast package in the experiment hub"))
      }else{
        keepIndVec = c(keepIndVec, keepInd) 
      }
    }
    breastData = breastData[keepIndVec]
  }

  
  
  esets <- list()
  #if()
  for(i in seq_len(length(breastData)))
  {
    dataName = breastData[i]$title
    esets[[i]] = breastData[[names(breastData)[i]]]
    names(esets)[i] = breastData[i]$title
  }
  
  ## -----------------------------------------------------------------------------
  ##Explicit removal of samples from specified datasets:
  ## -----------------------------------------------------------------------------
  delim <- ":"   ##This is the delimiter used to specify dataset:sample,
  
  ## same as used in metagx getbrcadata
  #load("inst\\extdata\\BenDuplicate.rda")
  #source(system.file("extdata", "patientselection.config", package="MetaGxBreast"))
  load(system.file("extdata", "duplicates.rda", package="MetaGxBreast"))
  
  rmix <- duplicates
  ii <- 1
  while (length(rmix) > ii){
    rmix <- rmix [!is.element(names(rmix), rmix[[ii]])]
    ii <- ii+1
  }
  rmix <- unique(unlist(rmix))
  
  message("Clean up the esets.")
    
  for (i in seq_len(length(esets))){
    eset <- esets[[i]]
    
    ##filter genes with standard deviation below the required quantile
    if(quantileCutoff > 0 && quantileCutoff < 1){
      eset <- filterQuantile(eset, q=quantileCutoff)
    }
    ##rescale to z-scores
    if(rescale == TRUE){
      Biobase::exprs(eset) <- t(scale(t(Biobase::exprs(eset))))
    }
    
    if(removeDuplicates == TRUE){
      keepix <- setdiff(colnames(eset@assayData$exprs), rmix)
      if(length(keepix) != length(colnames(eset@assayData$exprs)))
      {
        newEset = ExpressionSet(Biobase::exprs(eset)[, keepix, drop=FALSE])
        newEset@experimentData = eset@experimentData
        newEset@phenoData = eset@phenoData
        newEset@phenoData@data = Biobase::pData(eset)[keepix, , drop=FALSE]
        newEset@featureData = eset@featureData
        eset = newEset
      }
      #Biobase::exprs(eset) <- Biobase::exprs(eset)[, keepix, drop=FALSE]
      #Biobase::pData(eset) <- Biobase::pData(eset)[keepix, , drop=FALSE]
      
    }
    
    ##include study if it has enough samples and events:
    if (!is.na(minNumberEvents)
        && exists("minSampleSize") && !is.na(minSampleSize)
        && minNumberEvents > 0
        && sum(eset$vital_status == "deceased") < minNumberEvents
        || ncol(eset) < minSampleSize)
    {
      message(paste("excluding",
                    "(minNumberEvents or minSampleSize)"))
      next
    }
    if(nrow(eset) < minNumberGenes) {
      message(paste("excluding experiment hub dataset",breastData[i]$title,"(minNumberGenes)"))
      next
    }
    if(removeRetracted && length(grep("retracted", Biobase::experimentData(eset)@other$warnings$warnings)) > 0){
      message(paste("excluding experiment hub dataset",breastData[i]$title,"(removeRetracted)"))
      next
    }
    if(removeSubsets && length(grep("subset", Biobase::experimentData(eset)@other$warnings$warnings)) > 0){
      message(paste("excluding experiment hub dataset",breastData[i]$title,"(removeSubsets)"))
      next
    }
    message(paste("including experiment hub dataset",breastData[i]$title))
    ##    featureNames(eset) <- make.names(featureNames(eset))  ##should not do this, it is irreversible.
    esets[[i]] <- eset
    rm(eset)
  }
    
  
  ##optionally take the intersection of genes common to all platforms:
  if(keepCommonOnly){
    features.per.dataset <- lapply(esets, Biobase::featureNames)
    intersect.genes <- intersectMany(features.per.dataset)
    esets <- lapply(esets, function(eset){
      eset <- eset[intersect.genes, ]
      return(eset)
    })
  }
  
  ids.with.missing.data <- which(vapply(esets, function(X)
    sum(!complete.cases(Biobase::exprs(X))) > 0, numeric(1)) == 1)
  message(paste("Ids with missing data:", paste(names(ids.with.missing.data),
                                                collapse=", ")))
  
  if (length(ids.with.missing.data) > 0 && imputeMissing) {
    for (i in ids.with.missing.data) {
      Biobase::exprs(esets[[i]]) = impute::impute.knn(Biobase::exprs(esets[[i]]))$data
    }
  }
  
  retList = list(esets, duplicates)
  names(retList) = c("esets", "duplicates")
  return(retList)
}
