\name{ MAQC2 }
\alias{ MAQC2 }
\docType{data}
\title{ MAQC2 }
\description{ ExpressionSet for the  MAQC2  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=20064235 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194 
  PMIDs: 20064235 
  No abstract available.
  notes:
   summary:     
      It is possible to build multi-gene classifiers of clinical outcome. Predic
tion accuracy depends on training sample size and classification difficult
y. 
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
assayData: 20967 features, 230 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      230 character character 

alt_sample_name: 
   Length     Class      Mode 
      230 character character 

sample_type: 
tumor 
  230 

er: 
negative positive 
      89      141 

pgr: 
negative positive 
     126      104 

her2: 
negative positive 
     190       40 

N: 
  0   1 
 66 164 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   45.00   51.00   52.02   59.00   79.00 

grade: 
  1   2   3 
 13  94 123 

treatment: 
chemotherapy 
         230 

batch: 
MAQC2 
  230 

uncurated_author_metadata: 
   Length     Class      Mode 
      230 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194 }
\keyword{datasets}
