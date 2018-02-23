\name{ STK }
\alias{ STK }
\docType{data}
\title{ STK }
\description{ ExpressionSet for the  STK  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16280042 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse1456 
  PMIDs: 16280042 
  No abstract available.
  notes:
   summary:     
      Expression profiling was able to better predict prognosis compared to hist
ological staging. 
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... 244889_at (36178 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36178 features, 159 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      159 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0    67.0   136.0   138.3   208.5   277.0 

sample_type: 
tumor 
  159 

er: 
negative positive 
      29      130 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   31.0    48.0    56.0    57.8    68.5    87.0 

grade: 
   1    2    3 NA's 
  28   58   61   12 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     66    2022    2467    2234    2846    3099 

recurrence_status: 
norecurrence   recurrence 
         113           46 

treatment: 
  chemotherapy hormonotherapy      untreated 
            89             48             22 

batch: 
STK 
159 

uncurated_author_metadata: 
   Length     Class      Mode 
      159 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse1456 }
\keyword{datasets}
