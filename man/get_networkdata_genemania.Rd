% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-7_db_genemania.R
\name{get_networkdata_genemania}
\alias{get_networkdata_genemania}
\title{get_networkdata_genemania()}
\usage{
get_networkdata_genemania(
  species = "Homo_sapiens",
  version = "current",
  cache = TRUE,
  add_annotation = TRUE,
  ...
)
}
\arguments{
\item{species}{from which species does the data come from}

\item{version}{version of the data files in genemania}

\item{cache}{default value set to TRUE (automatically checks if the data file is already stored in the cache)}

\item{add_annotation}{expanding the dataframe with four columns (Entrez_ID and Ensembl_ID)}

\item{...}{further arguments passed to or from other methods}
}
\value{
ppis_genemania
}
\description{
get_networkdata_genemania()
}
\examples{
\donttest{
db_genemania_df <- get_networkdata_genemania(
  species = "Homo_sapiens",
  version = "current"
)

db_genemania_df
}
}
