\name{ LUND2 }
\alias{ LUND2 }
\docType{data}
\title{ LUND2 }
\description{ ExpressionSet for the  LUND2  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=17452630 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5325 
  PMIDs: 17452630 
  No abstract available.
  notes:
   summary:     
      Microarray signature was able to show PTEN mRNA losse_when IHC was unable,
 even though tumors exhibited PTEN loss behavior. Stathmim was an accurate
 IHC marker of the signature and had prognostic significance.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1 2 ... 27648 (22008 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 22008 features, 105 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      105 character character 

alt_sample_name: 
   Length     Class      Mode 
      105 character character 

sample_type: 
tumor 
  105 

er: 
negative positive 
      60       45 

treatment: 
hormonotherapy 
           105 

batch: 
LUND2 
  105 

uncurated_author_metadata: 
   Length     Class      Mode 
      105 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5325 }
\keyword{datasets}
