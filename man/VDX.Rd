\name{ VDX }
\alias{ VDX }
\docType{data}
\title{ VDX }
\description{ ExpressionSet for the  VDX  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information:  http://www.ncbi.nlm.nih.gov/pubmed/?term=15721472; http://www.ncbi.nlm.nih.gov/pubmed/?term=17420468 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5327 
  PMIDs: 15721472, 17420468 
  No abstract available.
  notes:
   summary:     
      15721472: A gene signature was identified that can accurately predict dist
ant metastases in node-negative cases. 17420468: Tumors with a lung metast
atic gene signature were shown to be larger.  
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (21169 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 21169 features, 344 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      344 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    3.0   122.8   605.5   575.7   836.5  2038.0 

sample_type: 
tumor 
  344 

er: 
negative positive 
     135      209 

tumor_size: 
   1    2    3    4 NA's 
 146  132    5    3   58 

N: 
  0 
344 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   44.00   52.00   53.88   63.00   83.00      58 

grade: 
   1    2    3 NA's 
   7   42  148  147 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     61    1254    2616    2377    3285    5201 

dmfs_status: 
norecurrence   recurrence 
         226          118 

treatment: 
untreated 
      344 

batch: 
 VDX VDXN 
 286   58 

uncurated_author_metadata: 
   Length     Class      Mode 
      344 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2034; http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse5327 }
\keyword{datasets}
