\name{ UPP }
\alias{ UPP }
\docType{data}
\title{ UPP }
\description{ ExpressionSet for the  UPP  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16141321 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3494 
  PMIDs: 16141321 
  No abstract available.
  notes:
   summary:     
      A 32-gene expression signature of p53 was identified that differentiates p
-53 mutant and wild-type tumors. The signature is more effective than sequ
ence-based assessments of p53 in predicting prognosis and therapeutic resp
onse.  
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
assayData: 36178 features, 251 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      251 character character 

alt_sample_name: 
   Length     Class      Mode 
      251 character character 

sample_type: 
tumor 
  251 

er: 
negative positive     NA's 
      34      213        4 

pgr: 
negative positive 
      61      190 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.200   1.500   2.000   2.243   2.562  13.000 

N: 
   0    1 NA's 
 158   84    9 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  28.00   52.00   64.00   62.11   72.00   93.00 

grade: 
   1    2    3 NA's 
  67  128   54    2 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     30    1870    3711    3007    3985    4654      17 

recurrence_status: 
norecurrence   recurrence         NA's 
         181           55           15 

treatment: 
hormonotherapy      untreated           NA's 
            80            142             29 

batch: 
UPPT UPPU 
  80  171 

uncurated_author_metadata: 
   Length     Class      Mode 
      251 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3494 }
\keyword{datasets}
