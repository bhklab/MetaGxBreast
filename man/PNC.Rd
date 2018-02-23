\name{ PNC }
\alias{ PNC }
\docType{data}
\title{ PNC }
\description{ ExpressionSet for the  PNC  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=21910250 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20711 
  PMIDs: GSE20711, PMID 21910250 
  No abstract available.
  notes:
   summary:     
      Breast tumors can be further divided than the currently known expression s
ubtypes based on DNA methylation profiles.
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
assayData: 42447 features, 92 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   4 observations deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
   88.0    25.0      NA    11.3      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       92 character character 

alt_sample_name: 
   Length     Class      Mode 
       92 character character 

sample_type: 
tumor 
   92 

er: 
negative positive     NA's 
      43       45        4 

pgr: 
negative positive     NA's 
      43       40        9 

her2: 
negative positive     NA's 
      64       26        2 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.900   1.700   2.500   2.758   3.000  10.000       6 

N: 
   0    1 NA's 
  43   40    9 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  32.16   48.57   53.90   55.97   64.84   82.13       4 

grade: 
   1    2    3 NA's 
  13    5   70    4 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     29     967    2216    2122    2931    5139       7 

recurrence_status: 
norecurrence   recurrence         NA's 
          49           36            7 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    318    1940    2372    2525    3043    5139       4 

vital_status: 
deceased   living     NA's 
      25       63        4 

batch: 
PNC 
 92 

uncurated_author_metadata: 
   Length     Class      Mode 
       92 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20711 }
\keyword{datasets}
