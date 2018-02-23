\name{ DUKE }
\alias{ DUKE }
\docType{data}
\title{ DUKE }
\description{ ExpressionSet for the  DUKE  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/16273092 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3143 
  PMIDs: 16273092 
  No abstract available.
  notes:
   summary:     
      It was shown that the activatione_status of several oncogenic pathways can
 bee_identified by gene expression signatures. These gene signatures ident
ify deregulation of pathways, associations with clinically relevant outcom
es, and characteristics of specific cancers and tumor subtypes.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1000_at 1001_at ... AFFX-MurIL4_at (12085 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12085 features, 171 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   1 observation deleted due to missingness 
      n  events  median 0.95LCL 0.95UCL 
 170.00   43.00    9.01    6.22      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      171 character character 

alt_sample_name: 
   Length     Class      Mode 
      171 character character 

sample_type: 
tumor 
  171 

er: 
negative positive 
      57      114 

pgr: 
negative positive     NA's 
      23       65       83 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.20    1.80    2.30    2.74    3.50    8.50      83 

N: 
   0    1 NA's 
  53   36   82 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  171.0   417.0   957.5  1235.0  1852.0  4069.0       1 

vital_status: 
deceased   living     NA's 
      43      127        1 

batch: 
DUKE 
 171 

uncurated_author_metadata: 
   Length     Class      Mode 
      171 character character 

duplicates: 
DUKE.DUKE_T00.622 DUKE.DUKE_T01.052 DUKE.DUKE_T01.522 DUKE.DUKE_T01.534 
                1                 1                 1                 1 
             NA's 
              167 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse3143 }
\keyword{datasets}
