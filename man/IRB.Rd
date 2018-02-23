\name{ IRB }
\alias{ IRB }
\docType{data}
\title{ IRB }
\description{ ExpressionSet for the  IRB  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/18297396 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5460 
  PMIDs: 18297396 
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
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (42447 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 42447 features, 129 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      129 character character 

alt_sample_name: 
   Length     Class      Mode 
      129 character character 

sample_type: 
tumor 
  129 

er: 
negative positive 
      53       76 

her2: 
negative positive 
      98       31 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.800   1.500   2.200   2.488   3.000   8.500 

N: 
 0  1 
64 65 

grade: 
 1  2  3 
27 32 70 

treatment: 
untreated 
      129 

batch: 
IRB 
129 

uncurated_author_metadata: 
   Length     Class      Mode 
      129 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5460 }
\keyword{datasets}
