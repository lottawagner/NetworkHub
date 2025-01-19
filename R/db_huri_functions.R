# get_networkdata_huri() -----------

#' get_networkdata_huri()
#'
#' @param species  from which species does the data come from (default human because currently only human data provided from huri)
#' @param type different datasets , more information on "http://www.interactome-atlas.org/about/" , c("HI-union", "Lit-BM")
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_huri
#' @return db_huri_ppi_anno_df
#' @return db_huri_anno_df
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_huri_df <- get_networkdata_huri(
#'   species = "human",
#'   type = "HI-union"
#' )
#'
#' db_huri_df
#' }

get_networkdata_huri <- function( species = "human",
                                  type = "HI-union",
                                  cache = TRUE,
                                  get_annotation = TRUE,
                                  add_annotation = TRUE,
                                  ...) {



  # list species is actualized for version huri "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in huri

  if (!(species %in% list_species_huri)) { # if species is not in the list
    stop("Species not found as specified by huri,",
         "please check some valid entries by running `list_species_huri`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "huri_",
    species,
    "_",
    type
  ) # definition of the resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" huri url
    huri_url <-
      urlmaker_huri(
        species = species,
        type = type
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_huri

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = huri_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_huri <- vroom::vroom(network_file, col_names = FALSE)
  #ppis_huri <- head(read.delim(network_file, sep = " "))

  #Ensembl
  colnames(ppis_huri)[colnames(ppis_huri) == "X1"] <- "Ensembl_A"
  colnames(ppis_huri)[colnames(ppis_huri) == "X2"] <- "Ensembl_B"

  #Interaction amount

  message ("HuRi currently provides this amount of interactions:")
  message(dim(ppis_huri))

  #add annotation
  annotation_db <-
    huri_db_annotations$anno_db_huri[match(species, huri_db_annotations$species)]

  if (get_annotation && !is.na(annotation_db)){

    db_huri_anno_df <- get_annotation_huri(ppi_huri = ppis_huri,
                                           species = species
                                           )

    message("...created annotation dataframe")

    if (add_annotation) {

      db_huri_ppi_anno_df <- add_annotation_huri(anno_df = db_huri_anno_df,
                                                 ppi_huri = ppis_huri,
                                                 species = species
      )
      message("...added annotation from *db_huri_anno_df* to *db_huri_ppi_anno_df*")

      return(db_huri_ppi_anno_df)
    }

    if (!add_annotation){
      return(db_huri_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_huri)
}


# outside of function ----------

list_species_huri <- c("human")

list_db_annotationdbi_huri <- c("org.Hs.eg.db")


huri_db_annotations <- data.frame(species = list_species_huri,
                                  anno_db_huri = list_db_annotationdbi_huri,
                                  row.names = list_species_huri
)

# get_annotation_huri() --------

#' get_annotation_huri ()
#'
#' @param species  from which species does the data come from
#' @param ppi_huri variable defined by ppis_huri in get_networkdata_huri()
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#'@return anno_df
#'@return ppi_huri
#'
#'@export
#'
#'
#' @examples
#' \dontrun{
#' db_huri_df <- get_networkdata_huri(species = "human",
#'                                    type = "HI-union",
#'                                    cache = TRUE,
#'                                    get_annotation = FALSE,
#'                                    add_annotation = FALSE
#'                                   )
#'
#' db_huri_anno_df <- get_annotation_huri(ppi_huri = db_huri_df,
#'                                        species = "human"
#'                                        )
#' }

get_annotation_huri <- function(ppi_huri,
                                species
                                ) {

  if (!(species %in% list_species_huri)) { # if species is not in the list
    stop("Species not found as specified by huri,",
         "please check some valid entries by running `list_species_huri`") # stop function and print
  }

  annotation_db <-
    huri_db_annotations$anno_db_huri[match(species, huri_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_huri$Ensembl_A, ppi_huri$Ensembl_B))
    anno_df <- data.frame(
      ensembl_id = all_prot_ids,
      genesymbol = AnnotationDbi::mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "ENSEMBL", column = "SYMBOL"),
      uniprot_id = AnnotationDbi::mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "ENSEMBL", column = "UNIPROT"),
      entrez_id = AnnotationDbi::mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "ENSEMBL", column = "ENTREZID"),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_huri(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_huri containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppi_huri)
  }
}



# add_annotation_huri() --------
#' add_annotation_huri ()
#'
#' @param anno_df annotation dataframe (for corresponding species in huri)
#' @param ppi_huri variable defined by ppis_huri in get_networkdata_huri()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_huri with annotation columns for each interactor (for corresponding species in huri)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_huri_df <- get_networkdata_huri(species = "human",
#'                                    type = "HI-union",
#'                                    cache = TRUE,
#'                                    get_annotation = FALSE,
#'                                    add_annotation = FALSE
#'                                   )
#'
#' db_huri_anno_df <- get_annotation_huri(ppi_huri = db_huri_df,
#'                                        species = "human"
#'                                        )
#'
#' db_huri_ppi_anno_df <- add_annotation_huri(ppi_huri = db_huri_df,
#'                                            anno_df = db_huri_anno_df,
#'                                            species = "human"
#'                                            )
#'
#'}

add_annotation_huri <- function(ppi_huri,
                                anno_df,
                                species) {

    #adding GeneSymbol
    ppi_huri$GeneSymbol_A <-
      anno_df$genesymbol[match(ppi_huri$Ensembl_A, anno_df$ensembl_id)]
    ppi_huri$GeneSymbol_B <-
      anno_df$genesymbol[match(ppi_huri$Ensembl_B, anno_df$ensembl_id)]

    #adding Uniprot
    ppi_huri$Uniprot_A <-
      anno_df$uniprot_id[match(ppi_huri$Ensembl_A, anno_df$ensembl_id)]
    ppi_huri$Uniprot_B <-
      anno_df$uniprot_id[match(ppi_huri$Ensembl_B, anno_df$ensembl_id)]

    #adding Entrez
    ppi_huri$Entrez_A <-
      anno_df$entrez_id[match(ppi_huri$Ensembl_A, anno_df$ensembl_id)]
    ppi_huri$Entrez_B <-
      anno_df$entrez_id[match(ppi_huri$Ensembl_B, anno_df$ensembl_id)]

    return(ppi_huri)
}


# build_graph_huri() -----

#' build_graph_huri()
#'
#' @param graph_data ppi data from huri
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
#' db_huri_df <- get_networkdata_huri(
#'   species = "human",
#'   type = "HI-union"
#' )
#'
#' db_huri_graph <- build_graph_huri(graph_data = db_huri_df,
#'                                         output_format = "igraph")
#' db_huri_graph #list of 9024
#' }
#'
#'


build_graph_huri <- function (graph_data,
                              output_format = "igraph",
                              min_score_threshold = NULL){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data_processed <- graph_data


  #check on dimension (amount of rows)
  dim(graph_data)
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


