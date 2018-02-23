\name{ FNCLCC }
\alias{ FNCLCC }
\docType{data}
\title{ FNCLCC }
\description{ ExpressionSet for the  FNCLCC  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=17659439 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7017 
  PMIDs: 17659439 
  No abstract available.
  notes:
   summary:     
      A potentially more powerful clinicogenomic model was created by combining 
a subset of relevant genes from an already published gene expression signa
ture and a commonly used clinical prognostic model (NPI). The genes in thi
s model are known to have a role in breast cancer, carcinogenesis, or chem
otherapy resistance.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: UMGC_00005 UMGC_00007 ... UMGC_09018 (6064 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 6064 features, 150 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      150 character character 

alt_sample_name: 
   Length     Class      Mode 
      150 character character 

sample_type: 
tumor 
  150 

N: 
  1 
150 

treatment: 
chemotherapy 
         150 

batch: 
FNCLCC 
   150 

uncurated_author_metadata: 
   Length     Class      Mode 
      150 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7017 }
\keyword{datasets}
