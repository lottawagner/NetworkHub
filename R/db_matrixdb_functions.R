# get_networkdata_matrixdb() -----------

#' get_networkdata_matrixdb()
#'
#' @param species  from which species does the data come from (default human because currently only human data provided from matrixdb)
#' @param type datasets provided by MatrixDB: "CORE" = MatrixDB manually curated interaction dataset
#' @param version version of the data files in MatrixDB
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_matrixdb
#' @return db_matrixdb_ppi_anno_df
#' @return db_matrixdb_anno_df
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_matrixdb_df <- get_networkdata_matrixdb(
#'   species = "human",
#'   type = "CORE",
#'   version = "4_0",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_matrixdb_df
#' }
get_networkdata_matrixdb <- function(species = "human",
                                     type = "CORE",
                                     version = "4_0",
                                     cache = TRUE,
                                     get_annotation = TRUE,
                                     add_annotation = TRUE,
                                     ...) {
  # MatrixDB 4.0 (2024-09-01)
  # UPDATEVERSION

  # check that the value for species is listed in matrixdb

  if (!(species %in% list_species_matrixdb)) { # if species is not in the list
    stop(
      "Species not found as specified by matrixdb,",
      "please check some valid entries by running `list_species_matrixdb`"
    ) # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "matrixdb_",
    species,
    "_",
    type,
    "_",
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
    # buildup from "base" matrixdb url
    matrixdb_url <-
      urlmaker_matrixdb(
        species = species,
        type = type,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_matrixdb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = matrixdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_matrixdb <- vroom::vroom(network_file, delim = "\t")
  # ppis_matrixdb <- head(read.delim(network_file, sep = ""))

  # Interaction amount
  message("matrixdb currently provides this amount of interactions:")
  message(dim(ppis_matrixdb))

  colnames(ppis_matrixdb)[colnames(ppis_matrixdb) == "#ID(s) interactor A"] <- "Identifier_A"
  colnames(ppis_matrixdb)[colnames(ppis_matrixdb) == "ID(s) interactor B"] <- "Identifier_B"

  # Adding Uniprot columns
  ppis_matrixdb$Uniprot_A <- ppis_matrixdb$Identifier_A
  ppis_matrixdb$Uniprot_B <- ppis_matrixdb$Identifier_B

  ppis_matrixdb$Uniprot_A <- ifelse(grepl("^uniprotkb:", ppis_matrixdb$Uniprot_A), ppis_matrixdb$Uniprot_A, NA)
  ppis_matrixdb$Uniprot_B <- ifelse(grepl("^uniprotkb:", ppis_matrixdb$Uniprot_B), ppis_matrixdb$Uniprot_B, NA)

  # UniProt
  ppis_matrixdb$Uniprot_A <- gsub("uniprotkb:", "", ppis_matrixdb$Uniprot_A)
  ppis_matrixdb$Uniprot_B <- gsub("uniprotkb:", "", ppis_matrixdb$Uniprot_B)

  # add annotation
  annotation_db <-
    matrixdb_db_annotations$anno_db_matrixdb[match(species, matrixdb_db_annotations$species)]

  if (get_annotation && !is.na(annotation_db)) {
    db_matrixdb_anno_df <- get_annotation_matrixdb(
      ppi_matrixdb = ppis_matrixdb,
      species = species
    )

    message("...created annotation dataframe")

    if (add_annotation) {
      db_matrixdb_ppi_anno_df <- add_annotation_matrixdb(
        anno_df = db_matrixdb_anno_df,
        ppi_matrixdb = ppis_matrixdb,
        species = species
      )

      message("...added annotation from *db_matrixdb_anno_df* to *db_matrixdb_ppi_anno_df*")

      return(db_matrixdb_ppi_anno_df)
    }

    if (!add_annotation) {
      return(db_matrixdb_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation) {
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_matrixdb)
}

# outside of function ----------

list_species_matrixdb <- c("human")

list_db_annotationdbi_matrixdb <- c("org.Hs.eg.db")


matrixdb_db_annotations <- data.frame(
  species = list_species_matrixdb,
  anno_db_matrixdb = list_db_annotationdbi_matrixdb,
  row.names = list_species_matrixdb
)

# get_annotation_matrixdb() --------

#' get_annotation_matrixdb ()
#'
#' @param ppi_matrixdb variable defined by ppis_matrixdb in get_networkdata_matrixdb()
#' @param species  from which species does the data come from
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @return ppi_matrixdb
#' @return anno_df
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_matrixdb_df <- get_networkdata_matrixdb(
#'   species = "human",
#'   type = "CORE",
#'   version = "4_0",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_matrixdb_anno_df <- get_annotation_matrixdb(
#'   ppi_matrixdb = db_matrixdb_df,
#'   species = "human"
#' )
#' }
get_annotation_matrixdb <- function(ppi_matrixdb,
                                    species) {
  if (!(species %in% list_species_matrixdb)) { # if species is not in the list
    stop(
      "Species not found as specified by matrixdb,",
      "please check some valid entries by running `list_species_matrixdb`"
    ) # stop function and print
  }

  annotation_db <-
    matrixdb_db_annotations$anno_db_matrixdb[match(species, matrixdb_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_matrixdb$Uniprot_A, ppi_matrixdb$Uniprot_B)[!is.na(c(ppi_matrixdb$Uniprot_A, ppi_matrixdb$Uniprot_B))])
    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      genesymbol = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL", na.rm = FALSE
      ),
      ensembl_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL", na.rm = FALSE
      ),
      entrez_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID", na.rm = FALSE
      ),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_matrixdb(..., add_annotation = FALSE, ...)\n",
      "You will get ppis_matrixdb containing annotation for Uniprot_ and GeneSymbol_."
    )
    return(ppi_matrixdb)
  }
}



# add_annotation_matrixdb() --------
#' add_annotation_matrixdb ()
#'
#' @param anno_df annotation dataframe (for corresponding species in matrixdb)
#' @param ppi_matrixdb variable defined by ppis_matrixdb in get_networkdata_matrixdb()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_matrixdb with annotation columns for each interactor (for corresponding species in matrixdb)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_matrixdb_df <- get_networkdata_matrixdb(
#'   species = "human",
#'   type = "CORE",
#'   version = "4_0",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_matrixdb_anno_df <- get_annotation_matrixdb(
#'   ppi_matrixdb = db_matrixdb_df,
#'   species = "human"
#' )
#'
#' db_matrixdb_ppi_anno_df <- add_annotation_matrixdb(
#'   ppi_matrixdb = db_matrixdb_df,
#'   anno_df = db_matrixdb_anno_df,
#'   species = "human"
#' )
#' }
add_annotation_matrixdb <- function(ppi_matrixdb,
                                    anno_df,
                                    species) {
  # adding GeneSymbol
  ppi_matrixdb$GeneSymbol_A <-
    anno_df$genesymbol[match(ppi_matrixdb$Uniprot_A, anno_df$uniprot_id)]
  ppi_matrixdb$GeneSymbol_B <-
    anno_df$genesymbol[match(ppi_matrixdb$Uniprot_B, anno_df$uniprot_id)]

  # adding Ensembl
  ppi_matrixdb$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_matrixdb$Uniprot_A, anno_df$uniprot_id)]
  ppi_matrixdb$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_matrixdb$Uniprot_B, anno_df$uniprot_id)]

  # adding Entrez
  ppi_matrixdb$Entrez_A <-
    anno_df$entrez_id[match(ppi_matrixdb$Uniprot_A, anno_df$uniprot_id)]
  ppi_matrixdb$Entrez_B <-
    anno_df$entrez_id[match(ppi_matrixdb$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_matrixdb)
}

# build_graph_matrixdb() -----

#' build_graph_matrixdb()
#'
#' @param graph_data ppi data from matrixdb
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
#' db_matrixdb_df <- get_networkdata_matrixdb(
#'   species = "human",
#'   type = "CORE"
#' )
#'
#' db_matrixdb_graph <- build_graph_matrixdb(
#'   graph_data = db_matrixdb_df,
#'   output_format = "igraph"
#' )
#' db_matrixdb_graph # list of 235
#' }
#'
build_graph_matrixdb <- function(graph_data,
                                 output_format = "igraph",
                                 min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data_processed <- graph_data

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
