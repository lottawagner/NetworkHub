% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-10_db_matrixdb.R
\name{annotation_matrixdb}
\alias{annotation_matrixdb}
\title{annotation_matrixdb ()}
\usage{
annotation_matrixdb(ppi_matrixdb, species, type, ...)
}
\arguments{
\item{ppi_matrixdb}{variable defined by ppis_matrixdb in get_networkdata_matrixdb()}

\item{species}{from which species does the data come from}

\item{type}{datasets provided by MatrixDB: "CORE" = MatrixDB manually curated interaction dataset}

\item{...}{further arguments passed to or from other methods}
}
\value{
ppis_matrixdb_annotated
}
\description{
annotation_matrixdb ()
}
\examples{
#\donttest{
#annotation_matrixdb(ppi_matrixdb, species = "human", type = "CORE")
#}
}
