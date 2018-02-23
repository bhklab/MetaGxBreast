\name{ TRANSBIG }
\alias{ TRANSBIG }
\docType{data}
\title{ TRANSBIG }
\description{ ExpressionSet for the  TRANSBIG  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=17545524 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gsE7390 
  PMIDs: 17545524 
  No abstract available.
  notes:
   summary:     
      The 76-gene signature was validated. The results supports the hypothesis t
hat utilizing the gene signature could reduce the number of patients who r
eceive unnecessary adjuvant therapy. 
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
assayData: 20967 features, 198 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

      n  events  median 0.95LCL 0.95UCL 
  198.0    56.0      NA    17.1      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      198 character character 

sample_type: 
tumor 
  198 

er: 
negative positive 
      64      134 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.600   1.700   2.000   2.181   2.500   5.000 

N: 
  0 
198 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  24.00   42.00   46.00   46.39   51.00   60.00 

grade: 
   1    2    3 NA's 
  30   83   83    2 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    125    2375    4384    3954    5566    9108 

dmfs_status: 
norecurrence   recurrence 
         147           51 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    121    1528    3534    3399    5130    8711 

recurrence_status: 
norecurrence   recurrence 
         112           86 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    146    2744    4562    4150    5610    9108 

vital_status: 
deceased   living 
      56      142 

treatment: 
untreated 
      198 

batch: 
VDXGUYU VDXIGRU  VDXKIU VDXOXFU  VDXRHU 
     36      50      51      24      37 

uncurated_author_metadata: 
   Length     Class      Mode 
      198 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gsE7390 }
\keyword{datasets}
