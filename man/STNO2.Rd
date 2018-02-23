\name{ STNO2 }
\alias{ STNO2 }
\docType{data}
\title{ STNO2 }
\description{ ExpressionSet for the  STNO2  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=12829800 
  Title:  
  URL: http://smd.princeton.edu/cgi-bin/publication/viewPublication.pl?pub_no=248 
  PMIDs: 12829800 
  No abstract available.
  notes:
   summary:     
      Distinct breast cancer subtypes were determined by gene expression profile
s and were validated in other published datasets.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: IMAGE:1020315 IMAGE:1030271 ... IMAGE:971399 (3663
    total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 3663 features, 118 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

      n  events  median 0.95LCL 0.95UCL 
 118.00   46.00    4.67    3.34      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      118 character character 

alt_sample_name: 
   Length     Class      Mode 
      118 character character 

sample_type: 
tumor 
  118 

er: 
negative positive     NA's 
      31       82        5 

tumor_size: 
   1    2    3    4 NA's 
   6   13   62   32    5 

N: 
   0    1 NA's 
  34   79    5 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   46.25   58.00   58.47   71.75   85.00 

grade: 
   1    2    3 NA's 
  11   49   53    5 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  122.0   396.0   761.0   927.9  1233.0  2800.0      23 

recurrence_status: 
norecurrence   recurrence 
          58           60 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     91     426     898    1019    1392    5722 

vital_status: 
deceased   living 
      46       72 

treatment: 
  chemotherapy hormonotherapy      untreated 
            23             73             22 

batch: 
STNO2 
  118 

uncurated_author_metadata: 
   Length     Class      Mode 
      118 character character 

}}
\source{ http://smd.princeton.edu/cgi-bin/publication/viewPublication.pl?pub_no=248 }
\keyword{datasets}
