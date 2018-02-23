\name{ UCSF }
\alias{ UCSF }
\docType{data}
\title{ UCSF }
\description{ ExpressionSet for the  UCSF  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=17428335; http://www.ncbi.nlm.nih.gov/pubmed/?term=14612510 
  Title:  
  URL: Not available 
  PMIDs: 17428335, 14612510  
  No abstract available.
  notes:
   summary:     
      A gene set was identified that correctly predicted outcomes more effective
ly than using histological markers.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: probe_1 probe_3 ... probe_10365 (8015 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 8015 features, 162 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   29 observations deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
 133.00   44.00   11.56    9.25      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      162 character character 

alt_sample_name: 
   Length     Class      Mode 
      162 character character 

sample_type: 
tumor 
  162 

er: 
negative positive     NA's 
      41      101       20 

pgr: 
negative positive     NA's 
      46       94       22 

her2: 
negative positive     NA's 
      35       19      108 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   1.800   2.000   2.682   3.200  11.000       7 

N: 
   0    1 NA's 
  67   82   13 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  28.00   44.00   53.00   56.61   70.00   88.00       9 

grade: 
   1    2    3 NA's 
  14   62   74   12 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47     897    2040    2084    2992    8267      29 

dmfs_status: 
norecurrence   recurrence 
         140           22 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47     861    1865    1985    2847    8267      29 

recurrence_status: 
norecurrence   recurrence 
         125           37 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47    1087    2054    2140    3087    8267      29 

vital_status: 
deceased   living     NA's 
      54       99        9 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated 
                31                 38                 61                 22 
              NA's 
                10 

batch: 
UCSF 
 162 

uncurated_author_metadata: 
   Length     Class      Mode 
      162 character character 

}}
\source{ Not available }
\keyword{datasets}
