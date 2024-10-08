% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-2_db_hint.R
\name{annotation_hint}
\alias{annotation_hint}
\title{annotation_hint ()}
\usage{
annotation_hint(ppi_hint, species, version, type)
}
\arguments{
\item{ppi_hint}{variable defined by ppis_hint in get_networkdata_hint()}

\item{species}{from which species does the data come from}

\item{version}{version of the data files in hint}

\item{type}{different interaction files provided by hint (all high-quality)}
}
\value{
ppis_hint_annotated
}
\description{
annotation_hint ()
}
\examples{
#\donttest{
# annotation_hint(ppi_hint, species = "HomoSapiens", version = "2024-06", type = "binary")
#}

}
