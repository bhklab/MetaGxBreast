\name{ METABRIC }
\alias{ METABRIC }
\docType{data}
\title{ METABRIC }
\description{ ExpressionSet for the  METABRIC  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/22522925 
  Title:  
  URL: https://www.ebi.ac.uk/ega/studies/EGAS00000000083 
  PMIDs: 22522925 
  No abstract available.
  notes:
   summary:     
      
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ILMN_1802380 ILMN_1736104 ... ILMN_1709472 (36155
    total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36155 features, 2136 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   165 observations deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
 1971.0   891.0    12.3    11.6    13.2 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
     2136 character character 

alt_sample_name: 
   Length     Class      Mode 
     2136 character character 

sample_type: 
healthy   tumor 
    144    1992 

er: 
negative positive     NA's 
     440     1508      188 

her2: 
negative positive     NA's 
     676      148     1312 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   1.700   2.300   2.621   3.000  18.200     164 

N: 
   0    1 NA's 
1042  950  144 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  21.93   51.36   61.78   61.13   70.76   96.29      13 

grade: 
   1    2    3 NA's 
 170  775  957  234 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      3    1498    2632    2948    4357    9218     147 

vital_status: 
deceased   living     NA's 
     891     1081      164 

treatment: 
chemo.plus.hormono       chemotherapy     hormonotherapy          untreated 
               196                226               1029                685 

batch: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  1.000   1.000   3.000   2.613   3.000   5.000     144 

uncurated_author_metadata: 
   Length     Class      Mode 
     2136 character character 

duplicates: 
   Length     Class      Mode 
     2136 character character 

}}
\source{ https://www.ebi.ac.uk/ega/studies/EGAS00000000083 }
\keyword{datasets}
