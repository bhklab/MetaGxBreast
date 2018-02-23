\name{ LUND }
\alias{ LUND }
\docType{data}
\title{ LUND }
\description{ ExpressionSet for the  LUND  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18430221 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31863 
  PMIDs: 18430221 
  No abstract available.
  notes:
   summary:     
      A significant difference was found between the ER positive subgroup and ER
 negative subgroup in the gene expression profiles.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: H200006618 H200006808 ... H300022925 (11154 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 11154 features, 143 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      143 character character 

alt_sample_name: 
   Length     Class      Mode 
      143 character character 

sample_type: 
tumor 
  143 

er: 
negative positive 
      29      114 

pgr: 
negative positive     NA's 
      47       88        8 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.200   1.100   1.500   1.486   1.800   4.000       2 

N: 
  0 
143 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  27.00   47.50   56.00   54.76   63.00   73.00 

batch: 
LUNDS1 LUNDS2 LUNDS3 LUNDS4 
    30     47     22     44 

uncurated_author_metadata: 
   Length     Class      Mode 
      143 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31863 }
\keyword{datasets}
