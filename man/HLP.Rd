\name{ HLP }
\alias{ HLP }
\docType{data}
\title{ HLP }
\description{ ExpressionSet for the  HLP  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=19688261 
  Title:  
  URL: http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-543/ 
  PMIDs: 19688261 
  No abstract available.
  notes:
   summary:     
      The results show evidence of different patterns of genetic aberrations in 
distinct molecular subtypes of breast cancer. Patterns of copy number aber
rations may drive biological phenomena characteristic to each subtype.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 9g8cQB1TZtuiix.ulU fJUdX0IAn_P9VLTgJU ...
    xopB7pPn18FJ067uDs (26536 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 26536 features, 53 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       53 character character 

alt_sample_name: 
   Length     Class      Mode 
       53 character character 

sample_type: 
tumor 
   53 

er: 
negative positive 
      28       25 

pgr: 
negative positive 
      33       20 

her2: 
negative positive 
      40       13 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  1.200   1.800   2.450   2.648   3.000   8.000       5 

N: 
   0    1 NA's 
  27   25    1 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  30.00   47.50   53.50   54.96   64.25   81.00       5 

grade: 
 3 
53 

batch: 
HLP 
 53 

uncurated_author_metadata: 
   Length     Class      Mode 
       53 character character 

}}
\source{ http://www.ebi.ac.uk/arrayexpress/experiments/E-TABM-543/ }
\keyword{datasets}
