# get_networkdata_mint() -----------

#' get_networkdata_mint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in mint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation default value set to TRUE (MINT annotation dataframe contains Uniprot, GeneSymbol and Ensembl)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_mint
#' @return ensembl_entry
#' @return ppi_mint_df_annotated
#'
#' @importFrom vroom vroom
#' @importFrom stringr str_extract
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_mint_df <- get_networkdata_mint(
#'   species = "Homo Sapiens",
#'   version = "current"
#' )
#'
#' db_mint_df
#' }
get_networkdata_mint <- function(species,
                                 version = "current",
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {
  # list species is actualized for version mint "current"
  # UPDATEVERSION

  # check that the value for species is listed in mint

  if (!(species %in% list_species_mint)) { # if species is not in the list
    stop(
      "Species not found as specified by mint,",
      "please check some valid entries by running `list_species_mint`"
    ) # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "mint_",
    species,
    "_v",
    version
  ) # definition of the resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" mint url
    mint_url <-
      urlmaker_mint(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_mint

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = mint_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_mint <- vroom::vroom(network_file, col_names = FALSE)
  # ppis_mint <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_mint))

  colnames(ppis_mint) <- c("Identifier_A", "Identifier_B", "Alternative identifier for interactor A", "Alternative identifier for interactor B", "Aliases for A", "Aliases for B", "method", "First author", "pubmed", "NCBI Taxonomy identifier for interactor A", "NCBI Taxonomy identifier for interactor B", "interaction_type", "Source databases and identifiers", "Interaction identifier(s) in the corresponding source database", "Confidence score")

  if (add_annotation) {
    ppi_mint_df_annotated <- ppis_mint

    # UniProt
    ppi_mint_df_annotated$Uniprot_A <- str_extract(ppi_mint_df_annotated$Identifier_A, "uniprotkb:([A-Z0-9]+)")
    ppi_mint_df_annotated$Uniprot_A <- gsub("uniprotkb:", "", ppi_mint_df_annotated$Identifier_A)

    ppi_mint_df_annotated$Uniprot_B <- str_extract(ppi_mint_df_annotated$Identifier_B, "uniprotkb:([A-Z0-9]+)")
    ppi_mint_df_annotated$Uniprot_B <- gsub("uniprotkb:", "", ppi_mint_df_annotated$Identifier_B)


    # GeneSymbol
    ppi_mint_df_annotated$GeneSymbol_A <- str_extract(ppi_mint_df_annotated$`Aliases for A`, "uniprotkb:([^\\(]+)\\(gene name\\)")
    ppi_mint_df_annotated$GeneSymbol_A <- gsub("uniprotkb:", "", ppi_mint_df_annotated$GeneSymbol_A)
    ppi_mint_df_annotated$GeneSymbol_A <- gsub("\\(gene name\\)", "", ppi_mint_df_annotated$GeneSymbol_A)

    ppi_mint_df_annotated$GeneSymbol_B <- str_extract(ppi_mint_df_annotated$`Aliases for B`, "uniprotkb:([^\\(]+)\\(gene name\\)")
    ppi_mint_df_annotated$GeneSymbol_B <- gsub("uniprotkb:", "", ppi_mint_df_annotated$GeneSymbol_B)
    ppi_mint_df_annotated$GeneSymbol_B <- gsub("\\(gene name\\)", "", ppi_mint_df_annotated$GeneSymbol_B)

    # Ensembl
    ppi_mint_df_annotated$Ensembl_A <- strsplit(ppi_mint_df_annotated$`Alternative identifier for interactor A`, "\\|") # split entry by "|" -> you get a list for each row
    ppi_mint_df_annotated$Ensembl_A <- lapply(ppi_mint_df_annotated$Ensembl_A, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
      ensembl_entry <- grep("^ensembl:", x, value = TRUE) # search for the entry (x) in every list that starts (^) with "ensembl"
      ensembl_entry <- gsub("ensembl:", "", ensembl_entry) # remove prefix "ensembl" by using gsub()
      return(ensembl_entry)
    })

    ppi_mint_df_annotated$Ensembl_B <- strsplit(ppi_mint_df_annotated$`Alternative identifier for interactor B`, "\\|") # split entry by "|" -> you get a list for each row
    ppi_mint_df_annotated$Ensembl_B <- lapply(ppi_mint_df_annotated$Ensembl_B, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
      ensembl_entry <- grep("^ensembl:", x, value = TRUE) # search for the entry (x) in every list that starts (^) with "ensembl"
      ensembl_entry <- gsub("ensembl:", "", ensembl_entry) # remove prefix "ensembl" by using gsub()
      return(ensembl_entry)
    })

    # Confidence score
    ppi_mint_df_annotated$`Confidence score` <- str_extract(ppi_mint_df_annotated$`Confidence score`, "intact-miscore:([0-9\\.]+)")
    ppi_mint_df_annotated$`Confidence score` <- gsub("intact-miscore:", "", ppi_mint_df_annotated$`Confidence score`)

    # Entrez - MINT doesn't provide info about entrez number

    return(ppi_mint_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_mint)
  }
}

# outside of function ----------

list_species_mint <- c(
  "Homo Sapiens",
  "Mus Musculus",
  "Drosophila Melanogaster",
  "Saccharomyces Cerevisiae"
)

# build_graph_mint() -----

#' build_graph_mint()
#'
#' @param graph_data ppi data from mint
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#'
#' @importFrom igraph graph.data.frame simplify
#'
#' @return my_graph
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_mint_df <- get_networkdata_mint(
#'   species = "Homo Sapiens",
#'   version = "current"
#' )
#'
#' db_mint_graph <- build_graph_mint(
#'   graph_data = db_mint_df,
#'   output_format = "igraph",
#'   min_score_threshold = "0.35"
#' )
#' db_mint_graph # list of 12010
#' }
#'
build_graph_mint <- function(graph_data,
                             output_format = "igraph",
                             min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$`Confidence score` <- as.numeric(graph_data$`Confidence score`)

  # Erstelle das Histogramm mit 50 bins (breaks)
  hist(graph_data$`Confidence score`, breaks = 50)

  # select ppi data >= minimal score
  if (!is.null(min_score_threshold)) {
    graph_data_processed <- graph_data[graph_data$`Confidence score` >= min_score_threshold, ]
  } else {
    graph_data_processed <- graph_data
  }

  # check on dimension (amount of rows)
  dim(graph_data)
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
