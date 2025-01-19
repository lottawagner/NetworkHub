# get_networkdata_biogrid() -----------

#' get_networkdata_biogrid()
#'
#' @param species  default value = "9606" - from which species does the data come from
#' @param version version of the data files in biogrid
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_biogrid
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_biogrid_df <- get_networkdata_biogrid(
#'   species = "9606",
#'   version = "4.4.238",
#'   cache = TRUE,
#'   add_annotation = FALSE
#' )
#' db_biogrid_df
#' }
get_networkdata_biogrid <- function(species = "9606",
                                    version = "4.4.238",
                                    cache = TRUE,
                                    add_annotation = TRUE,
                                    ...) {
  # list species is actualized for version biogrid "current" (sept 2024)
  # UPDATEVERSION

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "biogrid_",
    species,
    "_v_",
    version
  )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" biogrid url
    biogrid_url <-
      urlmaker_biogrid(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_biogrid

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = biogrid_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_biogrid <- vroom::vroom(network_file)
  # ppis_biogrid <- head(read.delim(network_file, sep = " "))

  # If you want to check all species contained in the file
  species_A <- unique(ppis_biogrid$`Organism ID Interactor A`)
  species_B <- unique(ppis_biogrid$`Organism ID Interactor B`)
  list_species_biogrid <- union(species_A, species_B)

  if (!(species %in% list_species_biogrid)) { # if species is not in the list
    stop(
      "Species not found as specified by biogrid,",
      "please check some valid entries by running `list_species_biogrid`"
    ) # stop function and print
  }

  # Filter dataframe for species
  ppis_biogrid_filtered <- ppis_biogrid[(ppis_biogrid$`Organism ID Interactor A` == species) &
    (ppis_biogrid$`Organism ID Interactor B` == species), ]

  # rename columns

  if (add_annotation) {
    # Entrez
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "Entrez Gene Interactor A"] <- "Entrez_A"
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "Entrez Gene Interactor B"] <- "Entrez_B"

    # Uniprot
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "SWISS-PROT Accessions Interactor A"] <- "Uniprot_A"
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "SWISS-PROT Accessions Interactor B"] <- "Uniprot_B"

    # GeneSymbol
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "Official Symbol Interactor A"] <- "GeneSymbol_A"
    colnames(ppis_biogrid_filtered)[colnames(ppis_biogrid_filtered) == "Official Symbol Interactor B"] <- "GeneSymbol_B"

    return(ppis_biogrid_filtered)
  }

  return(ppis_biogrid_filtered)
}

# build_graph_biogrid() -----

#' build_graph_biogrid()
#'
#' @param graph_data ppi data from biogrid
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#'
#' @importFrom igraph graph.data.frame simplify
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_biogrid_df <- get_networkdata_biogrid(
#'   species = "9606",
#'   version = "4.4.238"
#' )
#'
#' db_biogrid_graph <- build_graph_biogrid(
#'   graph_data = db_biogrid_df,
#'   output_format = "igraph"
#' )
#' db_biogrid_graph # list of 20124
#' }
#'
build_graph_biogrid <- function(graph_data,
                                output_format = "igraph",
                                min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data_processed <- graph_data

  # check on dimension (amount of rows)

  dim(graph_data_processed)

  edges <- data.frame(
    from = graph_data_processed$GeneSymbol_A,
    to = graph_data_processed$GeneSymbol_B
  )

  # Create unique nodes (combine both GeneSymbol columns)
  nodes <- data.frame(
    id = unique(c(
      graph_data_processed$GeneSymbol_A,
      graph_data_processed$GeneSymbol_B
    )),
    label = unique(c(
      graph_data_processed$GeneSymbol_A,
      graph_data_processed$GeneSymbol_B
    ))
  )

  # If output format is igraph, return the igraph object
  if (output_format == "igraph") {
    whole_graph <- igraph::graph.data.frame(d = edges, directed = FALSE)
    my_graph <- igraph::simplify(whole_graph)
    return(my_graph)
  }
  # simplify by avoiding multiple entries?
  ## could make it smaller and easier to handle, without losing too much/at all in info
}
