\name{ EORTC10994 }
\alias{ EORTC10994 }
\docType{data}
\title{ EORTC10994 }
\description{ ExpressionSet for the  EORTC10994  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=15897907 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1561 
  PMIDs: 15897907 
  No abstract available.
  notes:
   summary:     
      The tumors with an apocrine gene expression profile had strong histologica
l apocrine features. These tumors were androgen receptor positive and were
 all ER negative, creating further classifications of tumor cells based on
 steroid receptor activity- luminal which are ER and AR positive, basal th
at are ER and AR negative, and molecular apocrine that are ER negative and
 AR positive.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (20967 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20967 features, 49 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       49 character character 

alt_sample_name: 
   Length     Class      Mode 
       49 character character 

sample_type: 
tumor 
   49 

er: 
negative positive 
      22       27 

pgr: 
negative positive     NA's 
      29       18        2 

tumor_size: 
 1  2  3  4 
 4 23 14  8 

N: 
 0  1 
19 30 

grade: 
   1    2    3 NA's 
   4   22   20    3 

batch: 
EORTC10994 
        49 

uncurated_author_metadata: 
   Length     Class      Mode 
       49 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1561 }
\keyword{datasets}
