# get_networkdata_matrixdb() -----------

#' get_networkdata_matrixdb()
#'
#' @param species  from which species does the data come from (default human because currently only human data provided from matrixdb)
#' @param type datasets provided by MatrixDB: "CORE" = MatrixDB manually curated interaction dataset
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with (GeneSymbol, Uniprot ID and Entrez_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_matrixdb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_matrixdb_df <- get_networkdata_matrixdb(
#'   species = "human",
#'   type = "CORE"
#' )
#'
#' db_matrixdb_df
#'

get_networkdata_matrixdb <- function( species = "human",
                                      type = "CORE",
                                      cache = TRUE,
                                      add_annotation = TRUE,
                                      ...) {
  # MatrixDB 4.0 (2024-09-01)
  # UPDATEVERSION

  # check that the value for species is listed in matrixdb

  if (!(species %in% list_species_matrixdb)) { # if species is not in the list
    stop("Species not found as specified by matrixdb,",
         "please check some valid entries by running `list_species_matrixdb`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "matrixdb_",
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
    # buildup from "base" matrixdb url
    matrixdb_url <-
      urlmaker_matrixdb(
        species = species,
        type = type
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_matrixdb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = matrixdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_matrixdb <- vroom::vroom(network_file, delim = "\t")
  #ppis_matrixdb <- head(read.delim(network_file, sep = ""))

  #Interaction amount
  message ("matrixdb currently provides this amount of interactions:")
  message(dim(ppis_matrixdb))

  colnames(ppis_matrixdb)[colnames(ppis_matrixdb) == "#ID(s) interactor A"] <- "Identifier_A"
  colnames(ppis_matrixdb)[colnames(ppis_matrixdb) == "ID(s) interactor B"] <- "Identifier_B"

  #Adding Uniprot columns
  ppis_matrixdb$Uniprot_A <- ppis_matrixdb$Identifier_A
  ppis_matrixdb$Uniprot_B <- ppis_matrixdb$Identifier_B

  ppis_matrixdb$Uniprot_A <- ifelse(grepl("^uniprotkb:", ppis_matrixdb$Uniprot_A), ppis_matrixdb$Uniprot_A, NA)
  ppis_matrixdb$Uniprot_B <- ifelse(grepl("^uniprotkb:", ppis_matrixdb$Uniprot_B), ppis_matrixdb$Uniprot_B, NA)

  #UniProt
  ppis_matrixdb$Uniprot_A <- gsub("uniprotkb:", "", ppis_matrixdb$Uniprot_A)
  ppis_matrixdb$Uniprot_B <- gsub("uniprotkb:", "", ppis_matrixdb$Uniprot_B)

  #add annotation
  annotation_db <-
    matrixdb_db_annotations$anno_db_matrixdb[match(species, matrixdb_db_annotations$species)]

  if (add_annotation) {
    ppi_matrixdb_df_annotated <- annotation_matrixdb(ppi_matrixdb = ppis_matrixdb,
                                                     species = species,
                                                     type = type)

    return(ppi_matrixdb_df_annotated)

  }

  # no annotation
  if (!add_annotation) {
    return(ppis_matrixdb)
  }
}

# outside of function ----------

list_species_matrixdb <- c("human")

list_db_annotationdbi_matrixdb <- c("org.Hs.eg.db")


matrixdb_db_annotations <- data.frame(species = list_species_matrixdb,
                                 anno_db_matrixdb = list_db_annotationdbi_matrixdb,
                                 row.names = list_species_matrixdb
)

# annotation_matrixdb() --------

#' annotation_matrixdb ()
#'
#' @param ppi_matrixdb variable defined by ppis_matrixdb in get_networkdata_matrixdb()
#' @param species  from which species does the data come from
#' @param type datasets provided by MatrixDB: "CORE" = MatrixDB manually curated interaction dataset
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#'@return ppis_matrixdb_annotated
#'
#'@export
#'
#'
#' @examples
#'
#' # annotation_matrixdb(ppi_matrixdb, species = "human", type = "CORE")
#'

annotation_matrixdb <- function(ppi_matrixdb,
                                species,
                                type,
                                ...) {
  # find database on corresponding species

  if (!(species %in% list_species_matrixdb)) { # if species is not in the list
    stop("Species not found as specified by matrixdb,",
         "please check some valid entries by running `list_species_matrixdb`") # stop function and print
  }

  annotation_db <-
    matrixdb_db_annotations$anno_db_matrixdb[match(species, matrixdb_db_annotations$species)]


  all_prot_ids <- unique(c(ppi_matrixdb$Uniprot_A, ppi_matrixdb$Uniprot_B) [!is.na(c(ppi_matrixdb$Uniprot_A, ppi_matrixdb$Uniprot_B))])
  anno_df <- data.frame(
    uniprot_id = all_prot_ids,
    genesymbol = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL", na.rm = FALSE),
    ensembl_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL", na.rm = FALSE),
    entrez_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID", na.rm = FALSE),
    row.names = all_prot_ids
  )


  ppis_matrixdb_annotated <- ppi_matrixdb

  #adding GeneSymbol
  ppis_matrixdb_annotated$GeneSymbol_A <-
    anno_df$genesymbol[match(ppis_matrixdb_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_matrixdb_annotated$GeneSymbol_B <-
    anno_df$genesymbol[match(ppis_matrixdb_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Ensembl
  ppis_matrixdb_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_matrixdb_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_matrixdb_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_matrixdb_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppis_matrixdb_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_matrixdb_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_matrixdb_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_matrixdb_annotated$Uniprot_B, anno_df$uniprot_id)]

  return(ppis_matrixdb_annotated)

}










#
