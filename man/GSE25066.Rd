\name{ GSE25066 }
\alias{ GSE25066 }
\docType{data}
\title{ GSE25066 }
\description{ ExpressionSet for the  GSE25066  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information:  
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066 
  PMIDs: 21558518 
  No abstract available.
  notes:
   summary:     
      
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
assayData: 20967 features, 508 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      508 character character 

alt_sample_name: 
   Length     Class      Mode 
      508 character character 

sample_type: 
tumor 
  508 

er: 
negative positive     NA's 
     205      297        6 

pgr: 
negative positive     NA's 
     258      243        7 

her2: 
negative positive     NA's 
     485        6       17 

T: 
 T0  T1  T2  T3  T4 
  3  30 255 145  75 

N: 
  0   1 
157 351 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   24.0    42.0    49.0    49.8    58.0    75.0 

grade: 
   1    2    3    4 NA's 
  32  180  259   15   22 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    0.0   636.5   999.9  1088.0  1500.0  2717.0 

dmfs_status: 
norecurrence   recurrence 
         397          111 

batch: 
GSE25066 
     508 

uncurated_author_metadata: 
   Length     Class      Mode 
      508 character character 

chemosensitivity_prediction: 
Rx Insensitive   Rx Sensitive 
           339            169 

GGI_prediction: 
High  Low 
 336  172 

PAM50_prediction: 
 Basal   Her2   LumA   LumB Normal 
   189     37    160     78     44 

dlda30_prediction: 
pCR  RD 
196 312 

RCB_prediction: 
   RCB-0/I RCB-II/III 
       230        278 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066 }
\keyword{datasets}
