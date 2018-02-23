\name{ DFHCC2 }
\alias{ DFHCC2 }
\docType{data}
\title{ DFHCC2 }
\description{ Test the efficacy of treating TNBC with neoadjuvant cisplatin; explore biomarkers to identify predictors of response }
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2834466/ 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18864 
  PMIDs: 20100965 

  Abstract: A 16 word abstract is available. Use 'abstract' method.
  notes:
   summary:     
      A subset of the patients experienced a response induced by cisplatin and b
iomarkers were identified that could predict response to cisplatin.  
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... AFFX-HUMISGF3A/M97935_MB_at
    (42447 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 42447 features, 84 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       84 character character 

unique_patient_ID: 
   Length     Class      Mode 
       84 character character 

sample_type: 
tumor 
   84 

er: 
negative positive 
      53       31 

pgr: 
negative positive 
      53       31 

her2: 
negative positive 
      66       18 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  29.00   45.00   53.00   52.89   59.00   85.00 

grade: 
 1  2  3 
10 16 58 

treatment: 
chemotherapy 
          84 

batch: 
DFHCC2_CISPLATIN DFHCC2_REFERENCE 
              24               60 

uncurated_author_metadata: 
   Length     Class      Mode 
       84 character character 

duplicates: 
   Length     Class      Mode 
       84 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18864 }
\keyword{datasets}
