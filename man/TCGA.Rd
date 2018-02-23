\name{ TCGA }
\alias{ TCGA }
\docType{data}
\title{ TCGA }
\description{ ExpressionSet for the  TCGA  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/23000897 
  Title:  
  URL: http://cancergenome.nih.gov/ 
  PMIDs: 23000897 
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
  featureNames: A1BG A2M ... ARHGAP11A.2 (19504 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19504 features, 1073 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

      n  events  median 0.95LCL 0.95UCL 
1073.00  103.00   10.05    8.56   12.05 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
     1073 character character 

alt_sample_name: 
   Length     Class      Mode 
     1073 character character 

unique_patient_ID: 
   Length     Class      Mode 
     1073 character character 

sample_type: 
tumor 
 1073 

er: 
negative positive     NA's 
     233      790       50 

pgr: 
negative positive     NA's 
     334      686       53 

her2: 
negative positive     NA's 
     549      161      363 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   49.00   58.00   58.48   68.00   90.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   -7.0   137.0   412.0   820.3  1180.0  6796.0 

vital_status: 
deceased   living 
     103      970 

batch: 
TCGA 
1073 

uncurated_author_metadata: 
   Length     Class      Mode 
     1073 character character 

}}
\source{ http://cancergenome.nih.gov/ }
\keyword{datasets}
