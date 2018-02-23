\name{dupplicates}
\alias{dupplicates}
\docType{data}
\title{
a list containing the names of patients that are believed to be dulicates across
datasets
}
\description{
The object is a list where each element is a patient ID that is believed to be a
duplicate of a patient in another dataset. Patients are designated as duplicated if
they have Spearman correlations greater than or equal to 0.98 with other patient
expression profiles
}
\format{
  A list with 107 elements, each of which is a patient ID.
}
\keyword{datasets}
