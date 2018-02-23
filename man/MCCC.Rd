\name{ MCCC }
\alias{ MCCC }
\docType{data}
\title{ MCCC }
\description{ ExpressionSet for the  MCCC  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=19960244 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19177 
  PMIDs: 19960244 
  No abstract available.
  notes:
   summary:     
      Overall, expression and copy number profiling of familial tumors have show
n that the tumors show molecular heterogeneity similar to sporadic tumors 
and are defined by their molecular subtypes rather than BRCA1 or BRCA2 ger
mline mutation status. 
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: probe_10017 probe_10021 ... probe_7650767 (19048 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19048 features, 75 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       75 character character 

sample_type: 
tumor 
   75 

batch: 
MCCC 
  75 

uncurated_author_metadata: 
   Length     Class      Mode 
       75 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19177 }
\keyword{datasets}
