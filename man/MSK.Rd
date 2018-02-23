\name{ MSK }
\alias{ MSK }
\docType{data}
\title{ MSK }
\description{ ExpressionSet for the  MSK  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16049480 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2603 
  PMIDs: 16049480 
  No abstract available.
  notes:
   summary:     
      A set of genes were identified that mark and mediate metastasis to the lun
g. Some genes confer growth advantages to both the breast tumor and lung e
nvironment, while others contribute to aggressive growth specifically in t
he lung. 
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
assayData: 20967 features, 99 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       99 character character 

alt_sample_name: 
   Length     Class      Mode 
       99 character character 

sample_type: 
tumor 
   99 

er: 
negative positive 
      42       57 

pgr: 
negative positive     NA's 
      55       43        1 

her2: 
positive     NA's 
      85       14 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.100   2.450   3.200   3.624   4.300  10.000 

N: 
 0  1 
34 65 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  30.00   46.50   56.00   55.81   63.50   87.00 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    245    1279    1971    1888    2575    3924      17 

dmfs_status: 
norecurrence   recurrence         NA's 
          55           27           17 

batch: 
MSK 
 99 

uncurated_author_metadata: 
   Length     Class      Mode 
       99 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2603 }
\keyword{datasets}
