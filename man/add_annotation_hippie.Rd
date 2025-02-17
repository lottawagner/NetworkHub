% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_hippie_functions.R
\name{add_annotation_hippie}
\alias{add_annotation_hippie}
\title{add_annotation_hippie ()}
\usage{
add_annotation_hippie(ppi_hippie, anno_df, species)
}
\arguments{
\item{ppi_hippie}{variable defined by ppis_hippie in get_networkdata_hippie()}

\item{anno_df}{annotation dataframe (for corresponding species in hippie)}

\item{species}{from which species does the data come from}
}
\value{
ppi_hippie with annotation columns for each interactor (for corresponding species in hippie)
}
\description{
add_annotation_hippie ()
}
\examples{
\dontrun{

db_hippie_df <- get_networkdata_hippie(
  species = "Homo_sapiens",
  version = "current",
  cache = TRUE,
  get_annotation = FALSE,
  add_annotation = FALSE
)

db_hippie_anno_df <- get_annotation_hippie(
  species = "Homo sapiens",
  version = "current"
)

db_hippie_ppi_anno_df <- add_annotation_hippie(
  ppi_hippie = db_hippie_df,
  anno_df = db_hippie_anno_df,
  species = "Homo_sapiens"
)
}
}
