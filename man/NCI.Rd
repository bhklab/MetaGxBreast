\name{ NCI }
\alias{ NCI }
\docType{data}
\title{ NCI }
\description{ ExpressionSet for the  NCI  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=12917485 
  Title:  
  URL: Supplemental data from paper 
  PMIDs: 12917485 
  No abstract available.
  notes:
   summary:     
      Expression patterns were strongly associated with ER status, moderately as
sociated with grade, but not associated with menopausal state, node status
, or tumor size. Genes that were significantly associated with survival we
re identified.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: AF106966 AF217974 ... Y12473 (5154 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5154 features, 99 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
       99 character character 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21580   21610   21640   21650   21670   21830 

sample_type: 
tumor 
   99 

er: 
negative positive 
      34       65 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.80    1.80    2.50    2.82    3.00    8.00 

N: 
 0  1 
46 53 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  33.00   49.00   57.00   57.47   64.50   90.00 

grade: 
 1  2  3 
16 38 45 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      8     967    2057    1969    2930    4067 

recurrence_status: 
norecurrence   recurrence 
          54           45 

treatment: 
  chemotherapy hormonotherapy      untreated 
            10             78             11 

batch: 
NCI 
 99 

uncurated_author_metadata: 
   Length     Class      Mode 
       99 character character 

}}
\source{ Supplemental data from paper }
\keyword{datasets}
