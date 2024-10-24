% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03-1_db_stringdb.R
\name{get_networkdata_stringdb}
\alias{get_networkdata_stringdb}
\title{get_networkdata_stringdb()}
\usage{
get_networkdata_stringdb(
  species,
  version,
  cache = TRUE,
  add_annotation = TRUE,
  ...
)
}
\arguments{
\item{species}{from which species does the data come from}

\item{version}{version of the data files in stringdb}

\item{cache}{default value set to TRUE (automatically checks if the data file is already stored in the cache)}

\item{add_annotation}{default value set to TRUE}

\item{...}{further arguments passed to or from other methods}
}
\value{
ppis_stringdb
}
\description{
get_networkdata_stringdb()
}
\examples{
\dontrun{
db_string_df <- get_networkdata_stringdb(species = "Mus musculus",
                                         version = "12.0"
                                        )
db_string_df
}

}
