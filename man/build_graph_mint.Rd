% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_mint_functions.R
\name{build_graph_mint}
\alias{build_graph_mint}
\title{build_graph_mint()}
\usage{
build_graph_mint(
  graph_data,
  output_format = "igraph",
  min_score_threshold = NULL
)
}
\arguments{
\item{graph_data}{ppi data from mint}

\item{output_format}{selection of different graph functions that can be used}

\item{min_score_threshold}{select ppis that are "confident" depending on the scoretype/value}
}
\value{
my_graph
}
\description{
build_graph_mint()
}
\examples{
\dontrun{

db_mint_df <- get_networkdata_mint(
  species = "Homo Sapiens",
  version = "current"
)

db_mint_graph <- build_graph_mint(
  graph_data = db_mint_df,
  output_format = "igraph",
  min_score_threshold = "0.35"
)
db_mint_graph # list of 12010
}

}
