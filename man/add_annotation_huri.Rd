% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_huri_functions.R
\name{add_annotation_huri}
\alias{add_annotation_huri}
\title{add_annotation_huri ()}
\usage{
add_annotation_huri(ppi_huri, anno_df, species)
}
\arguments{
\item{ppi_huri}{variable defined by ppis_huri in get_networkdata_huri()}

\item{anno_df}{annotation dataframe (for corresponding species in huri)}

\item{species}{from which species does the data come from}
}
\value{
ppi_huri with annotation columns for each interactor (for corresponding species in huri)
}
\description{
add_annotation_huri ()
}
\examples{
\dontrun{
db_huri_df <- get_networkdata_huri(
  species = "human",
  type = "HI-union",
  cache = TRUE,
  get_annotation = FALSE,
  add_annotation = FALSE
)

db_huri_anno_df <- get_annotation_huri(
  ppi_huri = db_huri_df,
  species = "human"
)

db_huri_ppi_anno_df <- add_annotation_huri(
  ppi_huri = db_huri_df,
  anno_df = db_huri_anno_df,
  species = "human"
)
}
}
