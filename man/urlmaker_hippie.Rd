% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/urlmaker_functions.R
\name{urlmaker_hippie}
\alias{urlmaker_hippie}
\title{urlmaker_hippie()}
\usage{
urlmaker_hippie(species = "Homo_sapiens", version = "current")
}
\arguments{
\item{species}{default value = "Homo_sapiens", because this database only provides human data}

\item{version}{default value = "current", version of the database ... #UPDATEVERSION}
}
\value{
url returns the corresponding url set by params
}
\description{
urlmaker_hippie()
}
\examples{

url_hippie <- urlmaker_hippie(species = "Homo_sapiens",
                            version = "current")
url_hippie

}
