\name{ CAL }
\alias{ CAL }
\docType{data}
\title{ CAL }
\description{ ExpressionSet for the  CAL  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/17157792 
  Title:  
  URL: http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-158/ 
  PMIDs: 17157792 
  No abstract available.
  notes:
   summary:     
      Recurrent copy number abnormalities differ between tumor subtypes as defin
ed by gene expression patterns.  Accuracy of stratification by outcome can
 be improved by combining expression and copy number.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (21169 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 21169 features, 118 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   1 observation deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
 117.00   77.00    8.96    8.33    9.71 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      118 character character 

sample_type: 
tumor 
  118 

er: 
negative positive 
      43       75 

pgr: 
negative positive     NA's 
      51       66        1 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.300   1.675   2.300   2.729   3.500   7.500       2 

N: 
 0  1 
51 67 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  31.00   44.00   51.00   55.06   66.00   88.00       1 

grade: 
   1    2    3 NA's 
  10   42   61    5 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     767    2059    2094    3336    5183       1 

dmfs_status: 
norecurrence   recurrence         NA's 
          91           26            1 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      0     767    2059    2094    3336    5183       1 

recurrence_status: 
norecurrence   recurrence         NA's 
          81           36            1 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     47    1117    2234    2347    3504    5183       1 

vital_status: 
deceased   living 
      77       41 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated 
                25                 36                 40                 14 
              NA's 
                 3 

batch: 
CAL 
118 

uncurated_author_metadata: 
   Length     Class      Mode 
      118 character character 

}}
\source{ http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-158/ }
\keyword{datasets}
