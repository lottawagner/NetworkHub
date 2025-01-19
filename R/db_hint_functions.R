# get_networkdata_hint() -----------

#' get_networkdata_hint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in hint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param type different interaction files provided by hint (all high-quality)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_hint
#' @return db_hint_ppi_anno_df
#' @return db_hint_anno_df
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_hint_df <- get_networkdata_hint(
#'   species = "HomoSapiens",
#'   version = "2024-06",
#'   type = "binary",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#' }
get_networkdata_hint <- function(species,
                                 version,
                                 type = "binary", # c("binary", "cocomp", "lcb", "lcc")
                                 cache = TRUE,
                                 get_annotation = TRUE,
                                 add_annotation = TRUE,
                                 ...) {
  # list species is actualized for version HINT "2024-06"
  # UPDATEVERSION

  # check that the value for species is listed in HINT

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop(
      "Species not found as specified by HINT,",
      "please check some valid entries by running `list_species_hint`"
    ) # stop function and print
  }

  # type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "hint_",
    species,
    "_v_",
    version,
    "_type_",
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
    # buildup from "base" hint url
    hint_url <-
      urlmaker_hint(
        type = type,
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_hint

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = hint_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_hint <- vroom::vroom(network_file)
  # ppis_hint <- head(read.delim(network_file, sep = " "))

  message(dim(ppis_hint))

  colnames(ppis_hint)[colnames(ppis_hint) == "Gene_A"] <- "GeneSymbol_A"
  colnames(ppis_hint)[colnames(ppis_hint) == "Gene_B"] <- "GeneSymbol_B"

  if (get_annotation && is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_hint(..., add_annotation = FALSE, ...)\n",
      "You will get ppis_hint containing annotation for Uniprot_ and GeneSymbol_."
    )
    return(ppis_hint)
  }

  if (get_annotation && !is.na(annotation_db)) {
    db_hint_anno_df <- get_annotation_hint(
      ppi_hint = ppis_hint,
      species = species
    )

    message("...created annotation dataframe")

    if (add_annotation) {
      db_hint_ppi_anno_df <- add_annotation_hint(
        anno_df = db_hint_anno_df,
        ppi_hint = ppis_hint,
        species = species
      )
      message("...added annotation from *db_hint_anno_df* to *db_hint_ppi_anno_df*")

      return(db_hint_ppi_anno_df)
    }

    if (!add_annotation) {
      return(db_hint_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation) {
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_hint)
}

# outside of function ----------

list_species_hint <- c(
  "HomoSapiens",
  "SaccharomycesCerevisiae",
  "SchizosaccharomycesPombe",
  "MusMusculus",
  "DrosophilaMelanogaster",
  "CaenorhabditisElegans",
  "ArabidopsisThaliana",
  "EscherichiaColi",
  "RattusNorvegicus",
  "OryzaSativa"
)

list_db_annotationdbi_hint <- c(
  "org.Hs.eg.db",
  "org.Sc.sgd.db",
  NA,
  "org.Mm.eg.db",
  "org.Dm.eg.db",
  "org.Ce.eg.db",
  "org.At.tair.db",
  NA,
  "org.Rn.eg.db",
  NA
)


hint_db_annotations <- data.frame(
  species = list_species_hint,
  anno_db_hint = list_db_annotationdbi_hint,
  row.names = list_species_hint
)


# get_annotation_hint() --------

#' get_annotation_hint ()
#'
#' @param ppi_hint variable defined by ppis_hint in get_networkdata_hint()
#' @param species  from which species does the data come from
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @import org.Sc.sgd.db
#' @import org.Mm.eg.db
#' @import org.Dm.eg.db
#' @import org.Ce.eg.db
#' @import org.At.tair.db
#' @import org.Rn.eg.db
#'
#'
#' @return anno_df
#' @return ppi_hint
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_hint_df <- get_networkdata_hint(
#'   species = "HomoSapiens",
#'   version = "2024-06",
#'   type = "binary",
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_hint_anno_df <- get_annotation_hint(
#'   ppi_hint = db_hint_df,
#'   species = "HomoSapiens"
#' )
#' }
#'
get_annotation_hint <- function(ppi_hint,
                                species) {
  # find database on corresponding species

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop(
      "Species not found as specified by HINT,",
      "please check some valid entries by running `list_species_hint`"
    ) # stop function and print
  }

  annotation_db <-
    hint_db_annotations$anno_db_hint[match(species, hint_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_hint$Uniprot_A, ppi_hint$Uniprot_B))
    # all_gene_ids <- unique(c(ppi_hint$GeneSymbol_A, ppi_hint$GeneSymbol_B))

    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      # gene_symbol = all_gene_ids,
      ensembl_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"
      ),
      entrez_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"
      ),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_hint(..., add_annotation = FALSE, ...)\n",
      "You will get ppis_hint containing annotation for Uniprot_ and GeneSymbol_."
    )
    return(ppi_hint)
  }
}

# add_annotation_hint() --------
#' add_annotation_hint ()
#'
#' @param anno_df annotation dataframe (for corresponding species in hint)
#' @param ppi_hint variable defined by ppis_hint in get_networkdata_hint()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_hint with annotation columns for each interactor (for corresponding species in hint)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_hint_df <- get_networkdata_hint(
#'   species = "HomoSapiens",
#'   version = "2024-06",
#'   type = "binary",
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_hint_anno_df <- get_annotation_hint(
#'   ppi_hint = db_hint_df,
#'   species = "HomoSapiens"
#' )
#'
#' db_hint_ppi_anno_df <- add_annotation_hint(
#'   ppi_hint = db_hint_df,
#'   anno_df = db_hint_anno_df,
#'   species = "HomoSapiens"
#' )
#' }
add_annotation_hint <- function(ppi_hint,
                                anno_df,
                                species) {
  # adding Ensembl
  ppi_hint$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_hint$Uniprot_A, anno_df$uniprot_id)]
  ppi_hint$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_hint$Uniprot_B, anno_df$uniprot_id)]

  # adding Entrez
  ppi_hint$Entrez_A <-
    anno_df$entrez_id[match(ppi_hint$Uniprot_A, anno_df$uniprot_id)]
  ppi_hint$Entrez_B <-
    anno_df$entrez_id[match(ppi_hint$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_hint)
}


# build_graph_hint() -----

#' build_graph_hint()
#'
#' @param graph_data ppi data from HINT
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
#' db_hint_df <- get_networkdata_hint(
#'   species = "HomoSapiens",
#'   version = "2024-06",
#'   type = "binary",
#'   cache = TRUE,
#'   get_annotation = TRUE,
#'   add_annotation = TRUE
#' )
#'
#' db_hint_graph <- build_graph_hint(
#'   graph_data = db_hint_df,
#'   output_format = "igraph",
#'   min_score_threshold = NULL
#' )
#' db_hint_graph # list of ...TODO
#' }
#'
build_graph_hint <- function(graph_data,
                             output_format = "igraph",
                             min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)


  # select ppi data >= minimal score
  if (!is.null(min_score_threshold)) {
    graph_data_processed <- graph_data[graph_data$high_quality, ]
  }

  # select all ppi data
  else {
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
