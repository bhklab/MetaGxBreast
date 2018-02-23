\name{ UNC4 }
\alias{ UNC4 }
\docType{data}
\title{ UNC4 }
\description{ ExpressionSet for the  UNC4  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=20813035 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse18229 
  PMIDs: 20813035 
  No abstract available.
  notes:
   summary:     
      Clinically, this subtype is usually triple negative invasive ductal carcin
omas with a poor prognosis. Response to standard of care preoperative chem
otherapy is intermediate between basal-like and luminal tumors. The claudi
n-low subtype most closely resembles the mammary epithelial stem cell.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: probe.10 probe.12 ... probe.79701 (5420 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5420 features, 305 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   65 observations deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
 240.00   51.00    7.73    6.82      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      305 character character 

sample_type: 
tumor 
  305 

er: 
negative positive     NA's 
      99      154       52 

pgr: 
negative positive     NA's 
     126      109       70 

her2: 
negative positive     NA's 
     203       58       44 

tumor_size: 
   1  1.5    3    6 NA's 
  60    1  129   43   72 

N: 
   0    1 NA's 
 126  135   44 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  24.00   46.00   55.00   56.73   68.00   89.00      59 

grade: 
   1    2    3 NA's 
  25   80  138   62 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   30.0   450.0   750.0   954.3  1380.0  3540.0      64 

recurrence_status: 
norecurrence   recurrence         NA's 
         170           70           65 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     30     540     885    1104    1590    5190      65 

vital_status: 
deceased   living     NA's 
      51      189       65 

batch: 
UNC4 
 305 

uncurated_author_metadata: 
   Length     Class      Mode 
      305 character character 

duplicates: 
 KOO.KOO_KF_105 KOO.KOO_T01_514            NA's 
              1               1             303 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse18229 }
\keyword{datasets}
