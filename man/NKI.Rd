\name{ NKI }
\alias{ NKI }
\docType{data}
\title{ NKI }
\description{ ExpressionSet for the  NKI  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=12490681; http://www.ncbi.nlm.nih.gov/pubmed/?term=11823860 
  Title:  
  URL: Not available 
  PMIDs: 12490681, 11823860  
  No abstract available.
  notes:
   summary:     
      It was found that the gene expression profile that was studied was more po
werful in predicting outcome of disease in younger patients than using sta
ndard clinical and pathological criteria.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: Contig45645_RC Contig44916_RC ... Contig62037_RC (14960
    total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 14960 features, 337 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   42 observations deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
    295      79      NA      NA      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      337 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    4.0   123.0   215.0   214.1   312.0   404.0 

sample_type: 
tumor 
  337 

er: 
negative positive 
      88      249 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.200   1.500   2.000   2.241   2.800   5.500 

N: 
  0   1 
193 144 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   26.0    40.0    45.0    44.2    49.0    62.0 

grade: 
  1   2   3 
 79 109 149 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      9    1252    2414    2546    3602    6699      18 

dmfs_status: 
norecurrence   recurrence         NA's 
         210          109           18 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      9    1252    2414    2546    3602    6699      18 

recurrence_status: 
norecurrence   recurrence         NA's 
         210          109           18 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     20    1934    2637    2870    3763    6694      42 

vital_status: 
deceased   living     NA's 
      79      216       42 

treatment: 
  chemotherapy hormonotherapy      untreated 
            90             40            207 

batch: 
 NKI NKI2 
 117  220 

uncurated_author_metadata: 
   Length     Class      Mode 
      337 character character 

}}
\source{ Not available }
\keyword{datasets}
