# get_networkdata_irefindex() -----------

#' get_networkdata_irefindex()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in irefindex
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppi_irefindex
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_irefindex_df <- get_networkdata_irefindex(species = "Homo sapiens",
#'                                              version = "08-28-2023",
#'                                              get_annotation = FALSE,
#'                                              add_annotation = FALSE
#'                                              )
#' db_irefindex_df
#' }

get_networkdata_irefindex <- function(species = "Homo sapiens",
                                      version = "08-28-2023",
                                      cache = TRUE,
                                      get_annotation = TRUE,
                                      add_annotation = TRUE,
                                      ...) {



  # list species is actualized for version irefindex "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in irefindex

  if (!(species %in% list_species_irefindex)) { # if species is not in the list
    stop("Species not found as specified by irefindex,",
         "please check some valid entries by running `list_species_irefindex`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "irefindex_",
    species,
    "_v",
    version
  ) # definition of the resource name

  species_id <- irefindex_db_annotations[species, "species_id"]

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" irefindex url
    irefindex_url <-
      urlmaker_irefindex(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_irefindex

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = irefindex_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_irefindex <- vroom::vroom(network_file)
  #ppis_irefindex <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_irefindex))

  #Uniprot
  colnames(ppis_irefindex)[colnames(ppis_irefindex) == "#uidA"] <- "Uniprot_A"
  colnames(ppis_irefindex)[colnames(ppis_irefindex) == "uidB"] <- "Uniprot_B"
  #remove complex and refseq entries
  rows_to_remove_1 <- grep("^complex:", ppis_irefindex$Uniprot_A)
  rows_to_remove_2 <- grep("^refseq:", ppis_irefindex$Uniprot_A)
  ppis_irefindex <- ppis_irefindex[-rows_to_remove_1, ]
  ppis_irefindex <- ppis_irefindex[-rows_to_remove_2, ]

  #remove "uniprotdb:" description for each entry
  ppis_irefindex$Uniprot_A <- str_extract(ppis_irefindex$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
  ppis_irefindex$Uniprot_A <- gsub("uniprotkb:", "", ppis_irefindex$Uniprot_A)
  ppis_irefindex$Uniprot_B <- str_extract(ppis_irefindex$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
  ppis_irefindex$Uniprot_B <- gsub("uniprotkb:", "", ppis_irefindex$Uniprot_B)
  ppis_irefindex <- ppis_irefindex[!is.na(ppis_irefindex$Uniprot_A) & !is.na(ppis_irefindex$Uniprot_B), ]

  #Confidence score
  ppis_irefindex$confidence <- str_extract(ppis_irefindex$confidence, "lpr:([0-9\\.]+)")
  ppis_irefindex$confidence <- gsub("lpr:", "", ppis_irefindex$confidence)

  # match the annotation db with the corresponding species
  annotation_db <-
    irefindex_db_annotations$anno_db_irefindex[match(species, irefindex_db_annotations$species_irefindex)]

  if (get_annotation && is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_irefindex(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_irefindex containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppis_irefindex)
  }

  if (get_annotation && !is.na(annotation_db)){

    db_irefindex_anno_df <- get_annotation_irefindex(ppi_irefindex = ppis_irefindex,
                                                     species = species)

    message("...created annotation dataframe")

    if (add_annotation) {

      db_irefindex_ppi_anno_df <- add_annotation_irefindex(anno_df = db_irefindex_anno_df,
                                                           ppi_irefindex = ppis_irefindex,
                                                           species = species
                                                          )
      message("...added annotation from *db_irefindex_anno_df* to *db_irefindex_ppi_anno_df*")

      return(db_irefindex_ppi_anno_df)
    }

    if (!add_annotation){
      return(db_irefindex_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_irefindex)
}



# outside of function ----------

list_species_irefindex <- c ( "Homo sapiens",
                              "Mus musculus",
                              "Saccharomyces cerevisiae S288C",
                              "Escherichia",
                              "Rattus norvegicus",
                              "Saccharomyces cerevisiae",
                              "Drosophila melanogaster",
                              "Caenorhabditis elegans"
                            )

species_id_irefindex <- c(  "9606",
                            "10090",
                            "559292",
                            "562",
                            "10116",
                            "4932",
                            "7227",
                            "6239")




list_db_annotationdbi_irefindex <- c("org.Hs.eg.db",
                                     "org.Mm.eg.db",
                                     "org.Sc.sgd.db",
                                     "org.EcK12.eg.db",
                                     "org.Rn.eg.db",
                                     "org.Sc.sgd.db",
                                     "org.Dm.eg.db",
                                     "org.Ce.eg.db"
                                     )


irefindex_db_annotations <- data.frame(species_irefindex = list_species_irefindex,
                                       species_id = species_id_irefindex,
                                       anno_db_irefindex = list_db_annotationdbi_irefindex,
                                       row.names = list_species_irefindex
                                      )

# get_annotation_irefindex() --------

#' get_annotation_irefindex ()
#'
#' @param species  from which species does the data come from
#' @param ppi_irefindex variable defined by ppis_irefindex in get_networkdata_irefindex()
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Sc.sgd.db
#' @import org.EcK12.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#' @import org.Dm.eg.db
#' @import org.Ce.eg.db
#'
#'
#' @return ppi_irefindex
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_irefindex_df <- get_networkdata_irefindex(species = "Homo sapiens",
#'                                              version = "08-28-2023",
#'                                              get_annotation = FALSE,
#'                                              add_annotation = FALSE
#'                                              )
#'
#' db_irefindex_anno_df <- get_annotation_irefindex( ppi_irefindex = db_irefindex_df,
#'                                                   species = "Homo sapiens"
#'                                                 )
#' }



get_annotation_irefindex <- function(ppi_irefindex,
                                     species) {


  # find database for corresponding species
  annotation_db <- irefindex_db_annotations$anno_db_irefindex[match(species, irefindex_db_annotations$species_irefindex)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_irefindex$Uniprot_A, ppi_irefindex$Uniprot_B))

    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      genesymbol = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL"),
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
      entrez_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_irefindex(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_irefindex containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppi_irefindex)
  }
}

# add_annotation_irefindex() --------
#' add_annotation_irefindex ()
#'
#' @param anno_df annotation dataframe (for corresponding species in irefindex)
#' @param ppi_irefindex variable defined by ppis_irefindex in get_networkdata_irefindex()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_irefindex with annotation columns for each interactor (for corresponding species in irefindex)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_irefindex_df <- get_networkdata_irefindex(species = "Homo sapiens",
#'                                              version = "08-28-2023",
#'                                              get_annotation = FALSE,
#'                                              add_annotation = FALSE
#'                                              )
#'
#' db_irefindex_anno_df <- get_annotation_irefindex( ppi_irefindex = db_irefindex_df,
#'                                                   species = "Homo sapiens"
#'                                                 )
#'
#' db_irefindex_ppi_anno_df <- add_annotation_irefindex( ppi_irefindex = db_irefindex_df,
#'                                                       anno_df = db_irefindex_anno_df,
#'                                                       species = "Homo sapiens"
#'                                                     )
#'
#'}

add_annotation_irefindex <- function(ppi_irefindex,
                               anno_df,
                               species) {


  #adding GeneSymbol
  ppi_irefindex$GeneSymbol_A <-
    anno_df$genesymbol[match(ppi_irefindex$Uniprot_A, anno_df$uniprot_id)]
  ppi_irefindex$GeneSymbol_B <-
    anno_df$genesymbol[match(ppi_irefindex$Uniprot_B, anno_df$uniprot_id)]


  #adding Ensembl
  ppi_irefindex$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_irefindex$Uniprot_A, anno_df$uniprot_id)]
  ppi_irefindex$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_irefindex$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppi_irefindex$Entrez_A <-
    anno_df$entrez_id[match(ppi_irefindex$Uniprot_A, anno_df$uniprot_id)]
  ppi_irefindex$Entrez_B <-
    anno_df$entrez_id[match(ppi_irefindex$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_irefindex)


}



# build_graph_irefindex() -----

#' build_graph_irefindex()
#'
#' @param graph_data ppi data from irefindex
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident": lpr score (lowest PMID re-use)

#' @importFrom igraph graph.data.frame simplify
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_irefindex_df <- get_networkdata_irefindex(species = "Homo sapiens",
#'                                              version = "08-28-2023")
#'
#' db_irefindex_graph <- build_graph_irefindex(graph_data = db_irefindex_df,
#'                                             output_format = "igraph",
#'                                             min_score_threshold = "100")
#'
#' db_irefindex_graph #list of 5798 (score = 100)
#' }
#'
#'


build_graph_irefindex <- function (graph_data,
                                  output_format = "igraph",
                                  min_score_threshold = NULL){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$confidence <- as.numeric(graph_data$confidence)

  # Erstelle das Histogramm mit 50 bins (breaks)
  hist(graph_data$confidence, breaks = 50)

  #select ppi data >= minimal score
  if (!is.null(min_score_threshold)){
    graph_data_processed <- graph_data[graph_data$confidence <= min_score_threshold, ]
  } else {
    graph_data_processed <- graph_data
  }

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



