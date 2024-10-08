% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-9_db_huri.R
\name{annotation_huri}
\alias{annotation_huri}
\title{annotation_huri ()}
\usage{
annotation_huri(ppi_huri, species, type, add_annotation = TRUE, ...)
}
\arguments{
\item{ppi_huri}{variable defined by ppis_huri in get_networkdata_huri()}

\item{species}{from which species does the data come from}

\item{type}{different datasets , more information on "http://www.interactome-atlas.org/about/"}

\item{add_annotation}{expanding the dataframe with six columns (GeneSymbol, Uniprot ID and Entrez_ID)}

\item{...}{further arguments passed to or from other methods}
}
\value{
ppis_huri_annotated
}
\description{
annotation_huri ()
}
\examples{
# \donttest{
# annotation_huri(ppi_huri, species = "human", type = "HI-union")
#}
}
