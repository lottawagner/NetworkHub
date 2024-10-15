
#' ## info_species_stringdb() --------------------------------------------------------
#'
#' #' info_species_stringdb()
#' #'
#' #' @param version version of the data files in stringdb
#' #'
#' #' @return df_species returns the file containing info about the species annotation
#' #'
#' #' @importFrom utils read.delim
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' \dontrun{
#' #' info_species_stringdb()
#' #' }
#' info_species_stringdb <- function(version = "12.0"){
#'
#'   # use sprintf function to insert current version into the url
#'
#'   url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt",
#'                                   version)
#'
#'   # read.delim the data from the species text file (columns separated using a delimiter)
#'   df_species <- read.delim(url(url_species_stringdb))
#'
#'   # returns the datafile df_species
#'   return(df_species)
#' }

## get_accessoryinfo_stringdb()  -------------

#' Get accessory info from STRINGDB
#'
#' @param species from which species does the data come from
#' @param version version of the data files in stringdb
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return proteininfo_stringdb variable assigned to the datafile that contains info about protein annotation
#' @importFrom vroom vroom
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_stringdb_acc <- get_accessoryinfo_stringdb(species = "Homo sapiens",
#'                                        version = "12.0")
#' }

get_accessoryinfo_stringdb <- function(species = "Homo sapiens",
                                       version = "12.0",
                                       cache = TRUE,
                                       ...) {

  # species is matched to the id in the info file

  species_id <- NULL

  url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt",
                                  version)

  # read.delim the data from the species text file (columns separated using a delimiter)
  df_species <- read.delim(url(url_species_stringdb))

  info_species <- df_species # we find the information about species_name and species_id in the file info_species_stringdb
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
   match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]

  # rname of accessory info defined

  rname <- paste0(
    "stringdb_accessoryInfo_",
    species,
    "v_",
    version
  )

  # fetch (& cache) the accessory info_file

  ## fetch from the cache
  if (cache) {
    message("Trying to fetch from cache...")
    proteininfo_file <- fetch_NetworkHub(rname)
  }

  ## cache file and fetch from the cache

  if (!cache | is.null(proteininfo_file)) {
    message("Downloading to cache...")
    stringdb_url <- urlmaker_stringdb(
      type = "protein_info",
      species = species,
      version = version
    )

    proteininfo_file <- cache_NetworkHub(
      rname = rname,
      fpath = stringdb_url
    )

  }

  # reading in the fetched information file
  proteininfo_stringdb <- vroom::vroom(proteininfo_file)

  return(proteininfo_stringdb)
}

## create_annotation_from_stringdbaccessory() ----

#' create_annotation_from_stringdbaccessory()
#'
#' @param accessory_info #protein_info file from stringdb defined in get_networkhubdata_stringdb <- only use this function inside the other
#'
#' @return accessory_info_df dataframe that assignes each protein_id to ensemble or gene_symbol
#' @export
#'
#' @examples
#' #\dontrun{
#' db_stringdb_acc <- get_accessoryinfo_stringdb(species = "Homo sapiens",
#'                                        version = "12.0")
#'
#' create_annotation_from_stringdbaccessory(accessory_info = db_stringdb_acc)
#' }

create_annotation_from_stringdbaccessory <- function(accessory_info) {

  # create a dataframe that contains the unique string_protein_ids as protein_id and rname
  accessory_info_df <- data.frame(
    protein_id = unique(accessory_info$`#string_protein_id`),
    row.names = unique(accessory_info$`#string_protein_id`)
  )

  #rename columns

  accessory_info_df$protein_id <- str_extract(accessory_info_df$protein_id, "ENSP[0-9]+")
  #TODO fix ensp for all species

  return(accessory_info_df)
  }

## get_networkdata_stringdb() ---------

#' get_networkdata_stringdb()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in stringdb
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_stringdb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_string_df <- get_networkdata_stringdb(species = "Mus musculus",
#'                                          version = "12.0"
#'                                         )
#' db_string_df
#' }
#'
get_networkdata_stringdb <- function(species,
                                     version,
                                     cache = TRUE,
                                     add_annotation = TRUE,
                                     ...){

  # matching for the species...
  species_id <- NULL # looking for the corresponding id in the following

  url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt",
                                  version)

  # read.delim the data from the species text file (columns separated using a delimiter)
  df_species <- read.delim(url(url_species_stringdb))

  info_species <- df_species # we find the information about species_name and species_id in the file info_species_stringdb
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
    match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]


  if (is.null(species_id))
    stop("No match found!")

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "stringdb_",
    species,
    "_v",
    version
  ) # definition of the  resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" stringdb url
    stringdb_url <-
      urlmaker_stringdb(
        type = "PPI",
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_stringdb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = stringdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_stringdb <- vroom::vroom(network_file)
  ## ppis_stringdb <- read.delim(network_file, sep = " ")

  # id remap_identifiers = TRUE we use accessory functions to get ensemble id, gene_symbol and entrez
  if (add_annotation){
    accessory_info <- get_accessoryinfo_stringdb(
      species = species,
      version = version
    )
  # create a dataframe using the accessory functions
  accessory_info_df <-
    create_annotation_from_stringdbaccessory(accessory_info)



  #accessory_info_df$protein_id <- str_extract(accessory_info_df$protein_id, "ENSP[0-9]+")

  ppis_stringdb$ensemblid_1 <- accessory_info_df[ppis_stringdb$protein1, ][["ensembl_id"]]
  ppis_stringdb$ensemblid_2 <- accessory_info_df[ppis_stringdb$protein2, ][["ensembl_id"]]
  ppis_stringdb$gene_1 <- accessory_info_df[ppis_stringdb$protein1, ][["gene_symbol"]]
  ppis_stringdb$gene_2 <- accessory_info_df[ppis_stringdb$protein2, ][["gene_symbol"]]
  #ppis_stringdb$entrezid_1 <- accessory_info_df[ppis_stringdb$protein1, ][["entrez_id"]]
  #ppis_stringdb$entrezid_2 <- accessory_info_df[ppis_stringdb$protein2, ][["entrez_id"]]
  }

  # return ppis-stringdb containing remaped identifiers
  return(ppis_stringdb)
}



#' ## build_graph_stringdb() ----
#'
#' #' build_graph_stringdb()
#' #'
#' #' @param graph_data data describing the nodes and edges of the network
#' #' @param output_format vector that specifies possible output formats ("igraph", "grpahnel", "sif")
#' #' @param min_score_threshold optional threshold to filter edges based on their combination score
#' #' @param subset_nodes optional vector of nodes to extract only a part of the network
#' #'
#' #' @return my_graph object
#' #'
#' #' @importFrom igraph igraph
#' #' @export
#' #'
#' #' @examples
#' #'
#' #' network_all <- build_graph_stringdb(graph_data = graph_data,
#' #'                                     output_format = "igraph",
#' #'                                     min_score_threshold = 600)
#' #' TODO
#' build_graph_stringdb <- function(graph_data, # data describing the nodes and edges of the network
#'                                  output_format = c("igraph", "graphnel", "sif"), # vector that specifies possible output formats
#'                                  min_score_threshold = NULL, # optional threshold to filter edges based on their combination score
#'                                  subset_nodes = NULL) { # optional vector of nodes to extract only a part of the network
#'
#'   # graph data comes from the reading in function
#'   # now we need to check on the...
#'   colnames(graph_data) #returns the names of the coulmns of the corresponding dataframe
#'
#'
#'
#'   # create a histogram of the scores in the column combined_score in graph_data to define the threshold
#'   hist(graph_data$combined_score, breaks = 50)
#'   score_threshold <- 200
#'
#'
#'   # if there is a min_score_threshold #filter and only store edges that have a higher/same value as the threshold
#'   if (!is.null(min_score_threshold)){
#'     graph_data_processed <- graph_data[graph_data$combined_score >= min_score_threshold, ]
#'
#'   # no threshold ? # work with all edges
#'   } else {
#'     graph_data_processed <- graph_data
#'   }
#'
#'   dim(graph_data)
#'   dim(graph_data_processed)
#'
#'   # create an igraph object out of filtered data which is specified as undirected
#'   whole_graph <-
#'     igraph::graph.data.frame(d = graph_data_processed, directed = FALSE)
#'
#'
#'   # if there are subset nodes # create a subgraph which only contains edges and nodes from subset_nodes
#'   if (length(subset_nodes) > 0){
#'     my_graph <- igraph::induced_subgraph(
#'       whole_graph,
#'       which(V(whole_graph)$name %in% subset_nodes)
#'     )
#'
#'   # if there are no subset_nodes # create the whole graph
#'   } else {
#'     my_graph <- whole_graph
#'   }
#'
#'   # simplify graph by removing multiple edges and loops
#'   my_graph <- igraph::simplify(my_graph)
#'
#'   #return my_graph
#'   return(my_graph)
#'
#' }


## map experiment data to graph #TODO
## iplot_graph #TODO
