\name{ DUKE2 }
\alias{ DUKE2 }
\docType{data}
\title{ DUKE2 }
\description{ Predicting response with gene signature }
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/18024211 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse6861 
  PMIDs: 18024211 

  Abstract: A 5 word abstract is available. Use 'abstract' method.
  notes:
   summary:     
      Retraction in Lancet Feb 2011 (21277543); Regimen specific signatures were
 able to predict pathological complete response. Selecting patients with t
hese gene signataures could increase the proportion of patients with pCR t
han by basing clinical decisions on clinical factors. 
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1053_3p_at 117_3p_at ... X79510cds_3p_s_at (45490
    total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 45490 features, 160 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      160 character character 

alt_sample_name: 
   Length     Class      Mode 
      160 character character 

sample_type: 
tumor 
  160 

er: 
negative positive 
     123       37 

pgr: 
negative positive     NA's 
     133       25        2 

N: 
   0    1 NA's 
  58   95    7 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   43.00   49.00   49.41   56.00   70.00      35 

grade: 
   1    2    3 NA's 
   2   37   70   51 

treatment: 
chemotherapy 
         160 

batch: 
DUKE2 
  160 

uncurated_author_metadata: 
   Length     Class      Mode 
      160 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse6861 }
\keyword{datasets}
