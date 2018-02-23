\name{ MAINZ }
\alias{ MAINZ }
\docType{data}
\title{ MAINZ }
\description{ ExpressionSet for the  MAINZ  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18593943 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121 
  PMIDs: 18593943 
  No abstract available.
  notes:
   summary:     
      Poor prognosis is noted in tumors with low ER expression, showing the high
est level of proliferative activity. In some tumors with highly expressed 
B-cell or T-cell metagenes, metastases rarely occurred, even with high pro
liferation and low ER expression.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (20967 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20967 features, 200 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      200 character character 

alt_sample_name: 
   Length     Class      Mode 
      200 character character 

sample_type: 
tumor 
  200 

er: 
negative positive 
      38      162 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.100   1.500   2.000   2.070   2.425   6.000 

N: 
  0 
200 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  25.00   50.00   60.00   59.98   69.00   90.00 

grade: 
  1   2   3 
 29 136  35 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30    1905    2715    2816    3855    7200 

dmfs_status: 
norecurrence   recurrence 
         154           46 

treatment: 
untreated 
      200 

batch: 
MAINZ 
  200 

uncurated_author_metadata: 
   Length     Class      Mode 
      200 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121 }
\keyword{datasets}
