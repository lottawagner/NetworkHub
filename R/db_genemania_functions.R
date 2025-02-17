# get_networkdata_genemania() -----------

#' get_networkdata_genemania()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in genemania
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_genemania
#' @return db_genemania_ppi_anno_df
#' @return db_genemania_anno_df
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_genemania_df <- get_networkdata_genemania(
#'   species = "Homo_sapiens",
#'   version = "current",
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#' }
get_networkdata_genemania <- function(species = "Homo_sapiens",
                                      version = "current",
                                      cache = TRUE,
                                      get_annotation = TRUE,
                                      add_annotation = TRUE,
                                      ...) {
  # list species is actualized for version genemania "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in genemania

  if (!(species %in% list_species_genemania)) { # if species is not in the list
    stop(
      "Species not found as specified by genemania,",
      "please check some valid entries by running `list_species_genemania`"
    ) # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "genemania_",
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
    # buildup from "base" genemania url
    genemania_url <-
      urlmaker_genemania(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_genemania

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = genemania_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_genemania <- vroom::vroom(network_file)
  # ppis_genemania <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_genemania))

  # Uniprot
  colnames(ppis_genemania)[colnames(ppis_genemania) == "Gene_A"] <- "Uniprot_A"
  colnames(ppis_genemania)[colnames(ppis_genemania) == "Gene_B"] <- "Uniprot_B"

  annotation_db <-
    genemania_db_annotations$anno_db_genemania[match(species, genemania_db_annotations$species)]

  if (get_annotation && is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_genemania(..., add_annotation = FALSE, ...)\n"
    )
    return(ppis_genemania)
  }

  if (get_annotation && !is.na(annotation_db)) {
    db_genemania_anno_df <- get_annotation_genemania(
      ppi_genemania = ppis_genemania,
      species = species
    )

    message("...created annotation dataframe")

    if (add_annotation) {
      db_genemania_ppi_anno_df <- add_annotation_genemania(
        anno_df = db_genemania_anno_df,
        ppi_genemania = ppis_genemania,
        species = species
      )

      message("...added annotation from *db_genemania_anno_df* to *db_genemania_ppi_anno_df*")

      return(db_genemania_ppi_anno_df)
    }

    if (!add_annotation) {
      return(db_genemania_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation) {
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_genemania)
}



# outside of function ----------

list_species_genemania <- c(
  "Arabidopsis_thaliana",
  "Caenorhabditis_elegans",
  "Danio_rerio",
  "Drosophila_melanogaster",
  "Escherichia_coli",
  "Homo_sapiens",
  "Mus_musculus",
  "Rattus_norvegicus",
  "Saccharomyces_cerevisiae"
)

list_db_annotationdbi_genemania <- c(
  "org.At.tair.db",
  "org.Ce.eg.db",
  "org.Cf.eg.db",
  "org.Dm.eg.db",
  "org.EcK12.eg.db",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "org.Sc.sgd.db"
)


genemania_db_annotations <- data.frame(
  species = list_species_genemania,
  anno_db_genemania = list_db_annotationdbi_genemania,
  row.names = list_species_genemania
)


# get_annotation_genemania() --------

#' get_annotation_genemania ()
#'
#' @param ppi_genemania variable defined by ppis_genemania in get_networkdata_genemania()
#' @param species  from which species does the data come from
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.At.tair.db
#' @import org.Ce.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.EcK12.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#'
#' @return ppi_genemania
#' @return anno_df
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_genemania_df <- get_networkdata_genemania(
#'   species = "Homo_sapiens",
#'   version = "current",
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_genemania_anno_df <- get_annotation_genemania(
#'   ppi_genemania = db_genemania_df,
#'   species = "Homo_sapiens"
#' )
#' }
get_annotation_genemania <- function(ppi_genemania,
                                     species) {
  # find database on corresponding species

  if (!(species %in% list_species_genemania)) { # if species is not in the list
    stop(
      "Species not found as specified by genemania,",
      "please check some valid entries by running `list_species_genemania`"
    ) # stop function and print
  }

  annotation_db <-
    genemania_db_annotations$anno_db_genemania[match(species, genemania_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_genemania$Uniprot_A, ppi_genemania$Uniprot_B))
    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      gene_symbol = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL"
      ),
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
  }

  if (is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_genemania(..., add_annotation = FALSE, ...)\n",
      "You will get ppis_genemania containing annotation for Uniprot_ and GeneSymbol_."
    )
    return(ppi_genemania)
  }
}


# add_annotation_genemania() --------
#' add_annotation_genemania ()
#'
#' @param anno_df annotation dataframe (for corresponding species in genemania)
#' @param ppi_genemania variable defined by ppis_genemania in get_networkdata_genemania()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_genemania with annotation columns for each interactor (for corresponding species in genemania)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_genemania_df <- get_networkdata_genemania(
#'   species = "Homo_sapiens",
#'   version = "current",
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_genemania_anno_df <- get_annotation_genemania(
#'   ppi_genemania = db_genemania_df,
#'   species = "Homo_sapiens"
#' )
#'
#' db_genemania_ppi_anno_df <- add_annotation_genemania(
#'   ppi_genemania = db_genemania_df,
#'   anno_df = db_genemania_anno_df,
#'   species = "Homo_sapiens"
#' )
#' }
add_annotation_genemania <- function(ppi_genemania,
                                     anno_df,
                                     species) {
  # adding GeneSymbol
  ppi_genemania$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppi_genemania$Uniprot_A, anno_df$uniprot_id)]
  ppi_genemania$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppi_genemania$Uniprot_B, anno_df$uniprot_id)]

  # adding Ensembl
  ppi_genemania$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_genemania$Uniprot_A, anno_df$uniprot_id)]
  ppi_genemania$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_genemania$Uniprot_B, anno_df$uniprot_id)]

  # adding Entrez
  ppi_genemania$Entrez_A <-
    anno_df$entrez_id[match(ppi_genemania$Uniprot_A, anno_df$uniprot_id)]
  ppi_genemania$Entrez_B <-
    anno_df$entrez_id[match(ppi_genemania$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_genemania)
}

# output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B


# build_graph_genemania() -----

#' build_graph_genemania()
#'
#' @param graph_data ppi data from genemania
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
#' db_genemania_df <- get_networkdata_genemania(
#'   species = "Homo_sapiens",
#'   version = "current"
#' )
#'
#' db_genemania_graph <- build_graph_genemania(
#'   graph_data = db_genemania_df,
#'   output_format = "igraph",
#'   min_score_threshold = "0.00005"
#' )
#' db_genemania_graph # list of 18644
#' }
#'
build_graph_genemania <- function(graph_data,
                                  output_format = "igraph",
                                  min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$Weight <- as.numeric(graph_data$Weight)

  # Erstelle das Histogramm mit 50 bins (breaks)

  hist(graph_data$Weight, breaks = 5000, xlim = c(0, 0.0003))

  # select ppi data >= minimal score
  if (!is.null(min_score_threshold)) {
    graph_data_processed <- graph_data[graph_data$Weight >= min_score_threshold, ]
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
