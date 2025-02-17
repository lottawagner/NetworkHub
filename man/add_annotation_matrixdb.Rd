% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_matrixdb_functions.R
\name{add_annotation_matrixdb}
\alias{add_annotation_matrixdb}
\title{add_annotation_matrixdb ()}
\usage{
add_annotation_matrixdb(ppi_matrixdb, anno_df, species)
}
\arguments{
\item{ppi_matrixdb}{variable defined by ppis_matrixdb in get_networkdata_matrixdb()}

\item{anno_df}{annotation dataframe (for corresponding species in matrixdb)}

\item{species}{from which species does the data come from}
}
\value{
ppi_matrixdb with annotation columns for each interactor (for corresponding species in matrixdb)
}
\description{
add_annotation_matrixdb ()
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

db_matrixdb_ppi_anno_df <- add_annotation_matrixdb(
  ppi_matrixdb = db_matrixdb_df,
  anno_df = db_matrixdb_anno_df,
  species = "human"
)
}
}
