\name{ EXPO }
\alias{ EXPO }
\docType{data}
\title{ EXPO }
\description{ ExpressionSet for the  EXPO  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information:  
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109 
  PMIDs:  
  No abstract available.
  notes:
   summary:     
      N/A
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
assayData: 42447 features, 353 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      353 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1005   21640  101100  134700  215900  486200 

sample_type: 
tumor 
  353 

er: 
negative positive     NA's 
      85      161      107 

pgr: 
negative positive     NA's 
     114      129      110 

her2: 
negative positive     NA's 
     166       61      126 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  25.00   45.00   55.00   59.44   67.50   95.00       1 

grade: 
   1    2    3 NA's 
  32  114  151   56 

batch: 
EXPO 
 353 

uncurated_author_metadata: 
   Length     Class      Mode 
      353 character character 

duplicates: 
EXPO.EXPO_GSM53027 EXPO.EXPO_GSM53059               NA's 
                 1                  1                351 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109 }
\keyword{datasets}
