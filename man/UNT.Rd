\name{ UNT }
\alias{ UNT }
\docType{data}
\title{ UNT }
\description{ ExpressionSet for the  UNT  Dataset}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name:  
  Laboratory:  
  Contact information: http://www.ncbi.nlm.nih.gov/pubmed/?term=16478745; http://www.ncbi.nlm.nih.gov/pubmed/?term=17401012 
  Title:  
  URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2990 
  PMIDs: 16478745, 17401012  
  No abstract available.
  notes:
   summary:     
      A gene expression grading index (GGI) was developed. The GGI reclassified 
grade 2 patients into two groups with low and high risks of recurrence.
   mapping.method:     
      maxRowVariance
   mapping.group:     
      EntrezGene.ID
   preprocessing:     
      As published by original author.

featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: 1007_s_at 1053_at ... 244889_at (36084 total)
  varLabels: probeset gene EntrezGene.ID best_probe
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 36084 features, 133 samples
Platform type: 
--------------------------- 
Available sample meta-data: 
--------------------------- 

sample_name: 
   Length     Class      Mode 
      133 character character 

alt_sample_name: 
   Length     Class      Mode 
      133 character character 

sample_type: 
tumor 
  133 

er: 
negative positive     NA's 
      40       86        7 

pgr: 
negative positive     NA's 
       6       56       71 

tumor_size: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   1.200   1.900   1.892   2.300   6.000 

N: 
  0 
133 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  24.00   44.00   53.00   51.79   60.00   73.00 

grade: 
   1    2    3 NA's 
  32   51   29   21 

dmfs_days: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     61    1338    2809    2724    4078    5305 

dmfs_status: 
norecurrence   recurrence         NA's 
          97           28            8 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     61    1338    2675    2687    3912    5305 

recurrence_status: 
norecurrence   recurrence         NA's 
          76           49            8 

treatment: 
untreated 
      133 

batch: 
 KIU OXFU 
  64   69 

uncurated_author_metadata: 
   Length     Class      Mode 
      133 character character 

}}
\source{ http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse2990 }
\keyword{datasets}
