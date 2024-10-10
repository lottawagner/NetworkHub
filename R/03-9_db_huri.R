# get_networkdata_huri() -----------

#' get_networkdata_huri()
#'
#' @param species  from which species does the data come from (default human because currently only human data provided from huri)
#' @param type different datasets , more information on "http://www.interactome-atlas.org/about/" , c("HI-union", "Lit-BM")
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with six columns (GeneSymbol, Uniprot ID and Entrez_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_huri
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \donttest{
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

  if (add_annotation) {
    ppi_huri_df_annotated <- annotation_huri(ppi_huri = ppis_huri,
                                           species = species,
                                           type = type)
    return(ppi_huri_df_annotated)
  }

  # no annotation
  if (!add_annotation) {
    return(ppis_huri)
  }
}

# outside of function ----------

list_species_huri <- c("human")

list_db_annotationdbi_huri <- c("org.Hs.eg.db")


huri_db_annotations <- data.frame(species = list_species_huri,
                                 anno_db_huri = list_db_annotationdbi_huri,
                                 row.names = list_species_huri
)

# annotation_huri() --------

#' annotation_huri ()
#'
#' @param species  from which species does the data come from
#' @param type different datasets , more information on "http://www.interactome-atlas.org/about/"
#' @param ppi_huri variable defined by ppis_huri in get_networkdata_huri()
#' @param add_annotation expanding the dataframe with six columns (GeneSymbol, Uniprot ID and Entrez_ID)
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#'@return ppis_huri_annotated
#'
#'@export
#'
#'
#' @examples
#' # \donttest{
#' # annotation_huri(ppi_huri, species = "human", type = "HI-union")
#' #}

annotation_huri <- function(ppi_huri,
                           species,
                           type,
                           add_annotation = TRUE,
                           ...) {
  # find database on corresponding species

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


    ppis_huri_annotated <- ppi_huri

    #adding GeneSymbol
    ppis_huri_annotated$GeneSymbol_A <-
      anno_df$genesymbol[match(ppis_huri_annotated$Ensembl_A, anno_df$ensembl_id)]
    ppis_huri_annotated$GeneSymbol_B <-
      anno_df$genesymbol[match(ppis_huri_annotated$Ensembl_B, anno_df$ensembl_id)]

    #adding Uniprot
    ppis_huri_annotated$Uniprot_A <-
      anno_df$uniprot_id[match(ppis_huri_annotated$Ensembl_A, anno_df$ensembl_id)]
    ppis_huri_annotated$Uniprot_B <-
      anno_df$uniprot_id[match(ppis_huri_annotated$Ensembl_B, anno_df$ensembl_id)]

    #adding Entrez
    ppis_huri_annotated$Entrez_A <-
      anno_df$entrez_id[match(ppis_huri_annotated$Ensembl_A, anno_df$ensembl_id)]
    ppis_huri_annotated$Entrez_B <-
      anno_df$entrez_id[match(ppis_huri_annotated$Ensembl_B, anno_df$ensembl_id)]

    return(ppis_huri_annotated)
  }


  if (is.na(annotation_db)) {
    return(ppi_huri)
  }
}










#
