\name{ NCCS }
\alias{ NCCS }
\docType{data}
\title{ NCCS }
\description{ ExpressionSet for the  NCCS  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=18636107 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5364 
  PMIDs: 18636107 
  No abstract available.
  notes:
   summary:     
      48 genes were identified that displayed highly restricted levels of expres
sion in tumors compared to normal tissues. This was validated in 11 indepe
ndent cohorts of different cancer types.
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
assayData: 20967 features, 183 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      183 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    1.0    46.5    92.0    92.0   137.5   183.0 

sample_type: 
tumor 
  183 

batch: 
NCCS 
 183 

uncurated_author_metadata: 
   Length     Class      Mode 
      183 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5364 }
\keyword{datasets}
