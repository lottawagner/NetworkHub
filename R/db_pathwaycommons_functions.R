# get_networkdata_pathwaycommons() -----------

#' get_networkdata_pathwaycommons()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in pathwaycommons
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_pathwaycommons
#' @return db_pathwaycommons_ppi_anno_df
#' @return db_pathwaycommons_anno_df
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_pathwaycommons_df <- get_networkdata_pathwaycommons( species = "human",
#'                                                         version = "v12",
#'                                                         cache = TRUE,
#'                                                         get_annotation = FALSE,
#'                                                         add_annotation = FALSE
#'                                                         )
#' }

get_networkdata_pathwaycommons <- function( species = "human",
                                       version = "v12",
                                       cache = TRUE,
                                       get_annotation = TRUE,
                                       add_annotation = TRUE,
                                       ...) {



  # list species is actualized for version pathwaycommons "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in pathwaycommons

  if (!(species %in% list_species_pathwaycommons)) { # if species is not in the list
    stop("Species not found as specified by pathwaycommons,",
         "please check some valid entries by running `list_species_pathwaycommons`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "pathwaycommons_",
    species,
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
    # buildup from "base" pathwaycommons url
    pathwaycommons_url <-
      urlmaker_pathwaycommons(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_pathwaycommons

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = pathwaycommons_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_pathwaycommons <- vroom::vroom(network_file)
  #ppis_pathwaycommons <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_pathwaycommons))

  #GeneSymbol
  colnames(ppis_pathwaycommons)[colnames(ppis_pathwaycommons) == "PARTICIPANT_A"] <- "GeneSymbol_A"
  colnames(ppis_pathwaycommons)[colnames(ppis_pathwaycommons) == "PARTICIPANT_B"] <- "GeneSymbol_B"

  annotation_db <-
    pathwaycommons_db_annotations$anno_db_pathwaycommons[match(species, pathwaycommons_db_annotations$species)]

  if (get_annotation && !is.na(annotation_db)){

    db_pathwaycommons_anno_df <- get_annotation_pathwaycommons( ppi_pathwaycommons = ppis_pathwaycommons,
                                                                species = species
                                                              )

    message("...created annotation dataframe")

    if (add_annotation) {

      db_pathwaycommons_ppi_anno_df <- add_annotation_pathwaycommons( anno_df = db_pathwaycommons_anno_df,
                                                                      ppi_pathwaycommons = ppis_pathwaycommons,
                                                                      species = species
                                                                    )

      message("...added annotation from *db_pathwaycommons_anno_df* to *db_pathwaycommons_ppi_anno_df*")

      return(db_pathwaycommons_ppi_anno_df)
    }

    if (!add_annotation){
      return(db_pathwaycommons_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_pathwaycommons)
}

# outside of function ----------

list_species_pathwaycommons <- c("human")

list_db_annotationdbi_pathwaycommons <- c("org.Hs.eg.db")


pathwaycommons_db_annotations <- data.frame(species = list_species_pathwaycommons,
                                            anno_db_pathwaycommons = list_db_annotationdbi_pathwaycommons,
                                            row.names = list_species_pathwaycommons
                                            )


# get_annotation_pathwaycommons() --------

#' get_annotation_pathwaycommons ()
#'
#' @param ppi_pathwaycommons variable defined by ppis_pathwaycommons in get_networkdata_pathwaycommons()
#' @param species from which species does the data come from
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db

#'
#' @return ppi_pathwaycommons
#' @return anno_df
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_pathwaycommons_df <- get_networkdata_pathwaycommons( species = "human",
#'                                                         version = "v12",
#'                                                         cache = TRUE,
#'                                                         get_annotation = FALSE,
#'                                                         add_annotation = FALSE
#'                                                         )
#'
#' db_pathwaycommons_anno_df <- get_annotation_pathwaycommons( ppi_pathwaycommons = db_pathwaycommons_df,
#'                                                             species = "human")
#' }



get_annotation_pathwaycommons <- function(ppi_pathwaycommons,
                                          species) {
  # find database on corresponding species

  if (!(species %in% list_species_pathwaycommons)) { # if species is not in the list
    stop("Species not found as specified by pathwaycommons,",
         "please check some valid entries by running `list_species_pathwaycommons`") # stop function and print
  }

  annotation_db <-
    pathwaycommons_db_annotations$anno_db_pathwaycommons[match(species, pathwaycommons_db_annotations$species)]

  if (!is.na(annotation_db)) {

  all_prot_ids <- unique(c(ppi_pathwaycommons$GeneSymbol_A, ppi_pathwaycommons$GeneSymbol_B))
  anno_df <- data.frame(
    genesymbol = all_prot_ids,
    uniprot_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "UNIPROT"),
    ensembl_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "ENSEMBL"),
    entrez_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "ENTREZID"),
    row.names = all_prot_ids
    )

  return(anno_df)
  }

  if (is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_pathwaycommons(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_pathwaycommons containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppi_pathwaycommons)
  }
}


# add_annotation_pathwaycommons() --------
#' add_annotation_pathwaycommons ()
#'
#' @param anno_df annotation dataframe (for corresponding species in pathwaycommons)
#' @param ppi_pathwaycommons variable defined by ppis_pathwaycommons in get_networkdata_pathwaycommons()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_pathwaycommons with annotation columns for each interactor (for corresponding species in pathwaycommons)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_pathwaycommons_df <- get_networkdata_pathwaycommons( species = "human",
#'                                                         version = "v12",
#'                                                         cache = TRUE,
#'                                                         get_annotation = FALSE,
#'                                                         add_annotation = FALSE
#'                                                         )
#'
#' db_pathwaycommons_anno_df <- get_annotation_pathwaycommons( ppi_pathwaycommons = db_pathwaycommons_df,
#'                                                             species = "human")
#'
#' db_pathwaycommons_ppi_anno_df <- add_annotation_pathwaycommons( ppi_pathwaycommons = db_pathwaycommons_df,
#'                                                                 anno_df = db_pathwaycommons_anno_df,
#'                                                                 species = "human"
#'                                                               )
#'
#'}

add_annotation_pathwaycommons <- function(ppi_pathwaycommons,
                                    anno_df,
                                    species) {
  #adding Uniprot
  ppi_pathwaycommons$Uniprot_A <-
    anno_df$uniprot_id[match(ppi_pathwaycommons$GeneSymbol_A, anno_df$genesymbol)]
  ppi_pathwaycommons$Uniprot_B <-
    anno_df$uniprot_id[match(ppi_pathwaycommons$GeneSymbol_B, anno_df$genesymbol)]

  #adding Ensembl
  ppi_pathwaycommons$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_pathwaycommons$GeneSymbol_A, anno_df$genesymbol)]
  ppi_pathwaycommons$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_pathwaycommons$GeneSymbol_B, anno_df$genesymbol)]

  #adding Entrez
  ppi_pathwaycommons$Entrez_A <-
    anno_df$entrez_id[match(ppi_pathwaycommons$GeneSymbol_A, anno_df$genesymbol)]
  ppi_pathwaycommons$Entrez_B <-
    anno_df$entrez_id[match(ppi_pathwaycommons$GeneSymbol_B, anno_df$genesymbol)]

  return(ppi_pathwaycommons)

}










# build_graph_pathwaycommons() -----

#' build_graph_pathwaycommons()
#'
#' @param graph_data ppi data from pathwaycommons
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
#' db_pathwaycommons_df <- get_networkdata_pathwaycommons(
#'   species = "human",
#'   version = "v12"
#' )
#'
#' db_pathwaycommons_graph <- build_graph_pathwaycommons(graph_data = db_pathwaycommons_df,
#'                                         output_format = "igraph")
#' db_pathwaycommons_graph #list of 60434
#' }
#'
#'


build_graph_pathwaycommons <- function (graph_data,
                                        output_format = "igraph",
                                        min_score_threshold = NULL){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data_processed <- graph_data


  #check on dimension (amount of rows)
  dim(graph_data_processed)

  edges <- data.frame(from = graph_data_processed$GeneSymbol_A,
                      to = graph_data_processed$GeneSymbol_B)

  # Create unique nodes (combine both GeneSymbol columns)
  nodes <- data.frame(id = unique(c(graph_data_processed$GeneSymbol_A,
                                    graph_data_processed$GeneSymbol_B)),
                      label = unique(c(graph_data_processed$GeneSymbol_A,
                                       graph_data_processed$GeneSymbol_B)))

  # If output format is igraph, return the igraph object
  if (output_format == "igraph") {
    whole_graph <- igraph::graph.data.frame(d = edges, directed = FALSE)
    my_graph <- igraph::simplify(whole_graph)
    return(my_graph)
  }
  # simplify by avoiding multiple entries?
  ## could make it smaller and easier to handle, without losing too much/at all in info

}


