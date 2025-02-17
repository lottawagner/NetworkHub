% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_matrixdb_functions.R
\name{get_annotation_matrixdb}
\alias{get_annotation_matrixdb}
\title{get_annotation_matrixdb ()}
\usage{
get_annotation_matrixdb(ppi_matrixdb, species)
}
\arguments{
\item{ppi_matrixdb}{variable defined by ppis_matrixdb in get_networkdata_matrixdb()}

\item{species}{from which species does the data come from}
}
\value{
ppi_matrixdb

anno_df
}
\description{
get_annotation_matrixdb ()
}
\examples{
\dontrun{
db_matrixdb_df <- get_networkdata_matrixdb(
  species = "human",
  type = "CORE",
  version = "4_0",
  cache = TRUE,
  get_annotation = FALSE,
  add_annotation = FALSE
)

db_matrixdb_anno_df <- get_annotation_matrixdb(
  ppi_matrixdb = db_matrixdb_df,
  species = "human"
)
}
}
