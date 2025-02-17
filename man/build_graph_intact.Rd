% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_intact_functions.R
\name{build_graph_intact}
\alias{build_graph_intact}
\title{build_graph_intact()}
\usage{
build_graph_intact(
  graph_data,
  output_format = "igraph",
  min_score_threshold = NULL
)
}
\arguments{
\item{graph_data}{ppi data from intact}

\item{output_format}{selection of different graph functions that can be used}

\item{min_score_threshold}{select ppis that are "confident" depending on the scoretype/value}
}
\value{
my_graph
}
\description{
build_graph_intact()
}
\examples{
\dontrun{

db_intact_df <- get_networkdata_intact(
  species = "taxid:9606(human)|taxid:9606(Homo sapiens)",
  version = "current"
)

db_intact_graph <- build_graph_intact(
  graph_data = db_intact_df,
  output_format = "igraph",
  min_score_threshold = "0.35"
)
db_intact_graph # list of 17490
}

}
