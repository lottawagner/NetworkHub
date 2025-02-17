% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_hippie_functions.R
\name{get_networkdata_hippie}
\alias{get_networkdata_hippie}
\title{get_networkdata_hippie()}
\usage{
get_networkdata_hippie(
  species,
  version,
  cache = TRUE,
  get_annotation = TRUE,
  add_annotation = TRUE,
  ...
)
}
\arguments{
\item{species}{default value = "Homo_sapiens", because this database only provides human data}

\item{version}{default value = "current", version of the database ... #UPDATEVERSION}

\item{cache}{default value = TRUE, (automatically checks if the data file is already stored in the cache)}

\item{get_annotation}{default value = TRUE, creation of an annotation dataframe , default value set to TRUE}

\item{add_annotation}{default value = TRUE, adding annotation to ppi dataframe, default value set to TRUE}

\item{...}{further arguments passed to or from other methods}
}
\value{
ppis_hippie
}
\description{
get_networkdata_hippie()
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
db_hippie_df
}

}
