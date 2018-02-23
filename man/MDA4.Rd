\name{ MDA4 }
\alias{ MDA4 }
\docType{data}
\title{ MDA4 }
\description{ ExpressionSet for the  MDA4  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16896004 
  Title:  
  URL: http://bioinformatics.mdanderson.org/pubdata.html 
  PMIDs: 16896004 
  No abstract available.
  notes:
   summary:     
      The developed 30-probe set has high sensitivity and negative predictive va
lue, accurately identifying 12 out of 13 patients with pCR and 27 out of 2
8 patients with residual disease.
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
assayData: 21169 features, 129 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      129 character character 

unique_patient_ID: 
   Length     Class      Mode 
      129 character character 

sample_type: 
tumor 
  129 

er: 
negative positive     NA's 
      48       79        2 

pgr: 
negative positive     NA's 
      73       54        2 

her2: 
negative positive 
     114       15 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   0.500   1.800   2.162   3.000  10.000       8 

N: 
   0    1 NA's 
  59   62    8 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  28.00   43.00   51.00   51.43   61.00   73.00 

treatment: 
chemotherapy 
         129 

batch: 
MDA4 
 129 

uncurated_author_metadata: 
   Length     Class      Mode 
      129 character character 

duplicates: 
MDA4.MDA4_M207 MDA4.MDA4_M400           NA's 
             1              1            127 

}}
\source{ http://bioinformatics.mdanderson.org/pubdata.html }
\keyword{datasets}
