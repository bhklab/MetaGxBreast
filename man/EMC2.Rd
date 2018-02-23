\name{ EMC2 }
\alias{ EMC2 }
\docType{data}
\title{ EMC2 }
\description{ ExpressionSet for the  EMC2  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/19421193 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse12276 
  PMIDs: 19421193 
  No abstract available.
  notes:
   summary:     
      Genes were identified that may increase the ability of breast cancer cells
 to infiltrate the blood-brain barrier.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (42447 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 42447 features, 204 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      204 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00   51.75  102.50  102.50  153.20  204.00 

sample_type: 
tumor 
  204 

N: 
   0 NA's 
  48  156 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0     335     640     799    1098    3507 

dmfs_status: 
norecurrence   recurrence 
          19          185 

treatment: 
chemotherapy    untreated 
         156           48 

batch: 
EMC2 
 204 

uncurated_author_metadata: 
   Length     Class      Mode 
      204 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse12276 }
\keyword{datasets}
