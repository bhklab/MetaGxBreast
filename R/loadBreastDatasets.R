#' Function to load breast cancer SummarizedExperiment objectsfrom the Experiment Hub
#'
#' This function returns breast cancer datasets from the hub and a vector of patients from the datasets that are duplicates based on a spearman correlation > 0.98
#' @param rescale apply centering and scaling to the expression sets (default FALSE)
#' @param minNumberGenes an integer specifying to remove expression sets with less genes than this number (default 0)
#' @param minNumberEvents an integer specifying how man survival events must be in the dataset to keep the dataset (default 0)
#' @param minSampleSize an integer specifying the minimum number of patients required in a summarizedExperiment (default 0)
#' @param keepCommonOnly remove entrezIDs not common to all datasets (default FALSE)
#' @param imputeMissing remove patients from datasets with missing expression values
#' @param removeDuplicates remove patients with a Spearman correlation greater than or equal to 0.98 with other patient expression profiles (default TRUE)
#' @return a list with 2 elements. The First element named summarizedExperiments contains the datasets. The second element named duplicates contains
#' a vector with patient IDs for the duplicate patients (those with  Spearman correlation greater than or equal to 0.98 with other patient expression profiles).
#' @export
#' @importFrom Biobase esApply featureNames sampleNames exprs pData experimentData ExpressionSet
#' @importFrom lattice levelplot
#' @importFrom impute impute.knn
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom AnnotationHub query
#' @importFrom stats complete.cases sd quantile
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames colData rowData
#' @examples
#'
#' experimentsAndDups = loadBreastDatasets()


loadBreastDatasets = function(rescale = FALSE, minNumberGenes = 0,
                                minNumberEvents = 0, minSampleSize = 0,
                                keepCommonOnly = FALSE, imputeMissing = FALSE, removeDuplicates = FALSE)
{
  duplicates = NULL
  #if(getRversion() >= "2.15.1")  utils::globalVariables(c("."), add = F)
  ## -----------------------------------------------------------------------------
  ## needed functions
  ## -----------------------------------------------------------------------------

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


  ## -----------------------------------------------------------------------------
  ##load the summarizedExperiments
  ## -----------------------------------------------------------------------------

  hub = ExperimentHub::ExperimentHub()
  #AnnotationHub::possibleDates(hub)
  #query(eh, c("MetaGxOvarian", "SummarizedExperiment"))
  breastData = query(hub, c("MetaGxBreast", "SummarizedExperiment"))
  #pancreas issues: loading dataset 6, but missing /v1/ for 3 datasets
  dataList <- list()
  for(i in seq_len(length(breastData)))
  {
    dataList[[i]] = breastData[[names(breastData)[i]]]
    names(dataList)[i] = breastData[i]$title
  }
  names(dataList) = gsub("_sumexp", "", names(dataList), ignore.case = TRUE)

  ## -----------------------------------------------------------------------------
  ##Explicit removal of samples from specified datasets:
  ## -----------------------------------------------------------------------------
  delim <- ":"   ##This is the delimiter used to specify dataset:sample,

  ## same as used in metagx getbrcadata
  #load("inst\\extdata\\BenDuplicate.rda")
  #source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
  load(system.file("extdata", "duplicates.rda", package="MetaGxBreast"))

  rmix <- duplicates
  ii <- 1
  while (length(rmix) > ii){
    rmix <- rmix [!is.element(names(rmix), rmix[[ii]])]
    ii <- ii+1
  }
  rmix <- unique(unlist(rmix))
  rmix = substr(rmix, unlist(lapply(gregexpr("\\.", rmix), function(x) x[[1]][1]+1)), nchar(rmix))

  message("Clean up the summarizedExperiments")
  remInds = c()
  for (i in seq_len(length(dataList))){
    data <- dataList[[i]]
    include = TRUE
    ##rescale to z-scores
    if(rescale == TRUE){
      SummarizedExperiment::assay(data) = t(scale(t(SummarizedExperiment::assay(data))))
    }

    if(removeDuplicates == TRUE){
      keepix <- setdiff(colnames(SummarizedExperiment::assay(data)), rmix)
      if(length(keepix) != length(colnames(SummarizedExperiment::assay(data))))
      {
        keepix = which(!colnames(SummarizedExperiment::assay(data)) %in% rmix)
        data = data[ ,keepix]
      }
    }

    ##include study if it has enough samples and events:

    phenoData = colData(data)
    remData = TRUE
    if(nrow(phenoData) - sum(is.na(phenoData$vital_status)) >= minNumberEvents) remData = FALSE
    if(nrow(phenoData) - sum(is.na(phenoData$recurrence_status)) >= minNumberEvents) remData = FALSE
    if(nrow(phenoData) >= minSampleSize) remData = FALSE

    if(remData == TRUE){
      message(paste("excluding", names(dataList)[i], "(minNumberEvents or minSampleSize)"))
      remInds = c(remInds, i)
      include = FALSE
    }


    if(nrow(data) < minNumberGenes) {
      message(paste("excluding experiment hub dataset",names(dataList)[i],"(minNumberGenes)"))
      remInds = c(remInds, i)
      include = FALSE

    }

    if(imputeMissing == TRUE){
      notNaInds = which(colSums(is.na(SummarizedExperiment::assay(data))) == 0)
      data = data[, notNaInds]
      if(length(notNaInds) == 0){
        message(paste("excluding experiment hub dataset",names(dataList)[i],
                      "as every patient has at least 1 NA expression value (imputmissing = TRUE)"))
        remInds = c(remInds, i)
        include = FALSE

      }
    }
    if(include == TRUE) message(paste("including experiment hub dataset",names(dataList)[i]))
    ##    featureNames(eset) <- make.names(featureNames(eset))  ##should not do this, it is irreversible.
    dataList[[i]] <- data
    rm(data)
  }
  if(length(remInds) > 0) dataList[unique(remInds)] = NULL


  ##optionally take the intersection of genes common to all platforms:
  if(keepCommonOnly & length(dataList) > 0){
    features.per.dataset <- lapply(dataList, function(x) rowData(x)$EntrezGene.ID)
    intersect.genes <- intersectMany(features.per.dataset)
    dataList <- lapply(dataList, function(data){
      data <- data[which(rowData(data)$EntrezGene.ID %in% intersect.genes), ]
      return(data)
    })
  }
  if(length(dataList) == 0) warning("input values resulted in no datasets being returned")


  retList = list(dataList, duplicates)
  names(retList) = c("summarizedExperiments", "duplicates")
  return(retList)
}
