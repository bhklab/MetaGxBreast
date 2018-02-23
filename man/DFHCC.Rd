\name{ DFHCC }
\alias{ DFHCC }
\docType{data}
\title{ DFHCC }
\description{ ExpressionSet for the  DFHCC  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2826790/ 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19615 
  PMIDs: 20098429 
  No abstract available.
  notes:
   summary:     
      A small number of over-expressed and over-amplified genes were significant
ly associated with early recurrence despite adjuvant therapy.  This was ve
rified in independent cohorts. 
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
assayData: 42447 features, 115 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      115 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    6.0   155.0   230.0   293.3   398.5   828.0 

sample_type: 
tumor 
  115 

er: 
negative positive 
      45       70 

pgr: 
negative positive 
      51       64 

her2: 
negative positive 
      79       36 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.800   1.350   2.100   2.312   2.850   6.500 

N: 
 0  1 
62 53 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  32.00   45.00   53.00   53.89   60.00   85.00 

grade: 
 1  2  3 
23 28 64 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30    1500    1920    1799    2325    2640 

dmfs_status: 
norecurrence   recurrence 
         101           14 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated 
                42                 38                 22                  7 
              NA's 
                 6 

batch: 
DFHCC 
  115 

uncurated_author_metadata: 
   Length     Class      Mode 
      115 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19615 }
\keyword{datasets}
