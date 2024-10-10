
# get_networkdata_hint() -----------

#' get_networkdata_hint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in hint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param type different interaction files provided by hint (all high-quality)
#' @param add_annotation expanding the dataframe with four columns ( Entrez_ID and Ensembl_ID )
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_hint
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \donttest{
#' db_hint_df <- get_networkdata_hint(
#'   species = "HomoSapiens",
#'   version = "2024-06",
#'   type = "binary"
#' )
#'
#' db_hint_df
#' }

get_networkdata_hint <- function(species,
                                 version,
                                 type = "binary", #c("binary", "cocomp", "lcb", "lcc")
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {



  # list species is actualized for version HINT "2024-06"
  # UPDATEVERSION

  # check that the value for species is listed in HINT

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop("Species not found as specified by HINT,",
         "please check some valid entries by running `list_species_hint`") # stop function and print
  }


  # type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "hint_",
    species,
    "_v",
    version,
    "_type",
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
  #ppis_hint <- head(read.delim(network_file, sep = " "))
#
  message(dim(ppis_hint))

  colnames(ppis_hint)[colnames(ppis_hint) == "Gene_A"] <- "GeneSymbol_A"
  colnames(ppis_hint)[colnames(ppis_hint) == "Gene_B"] <- "GeneSymbol_B"

  annotation_db <-
    hint_db_annotations$anno_db_hint[match(species, hint_db_annotations$species)]

  if (add_annotation && !is.na(annotation_db)) {
    ppi_hint_df_annotated <- annotation_hint(ppi_hint = ppis_hint,
                                             species = species,
                                             version = version,
                                             type = type)
    return(ppi_hint_df_annotated)
  }

  if (add_annotation && is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_hint(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_hint containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppis_hint)
  }

  if (!add_annotation) {
    return(ppis_hint)
  }

}

# outside of function ----------

list_species_hint <- c("HomoSapiens",
                       "SaccharomycesCerevisiae",
                       "SchizosaccharomycesPombe",
                       "MusMusculus",
                       "DrosophilaMelanogaster",
                       "CaenorhabditisElegans",
                       "ArabidopsisThaliana",
                       "EscherichiaColi",
                       "RattusNorvegicus",
                       "OryzaSativa")

list_db_annotationdbi_hint <- c("org.Hs.eg.db",
                               "org.Sc.sgd.db",
                               NA,
                               "org.Mm.eg.db",
                               "org.Dm.eg.db",
                               "org.Ce.eg.db",
                               "org.At.tair.db",
                               NA,
                               "org.Rn.eg.db",
                               NA)


hint_db_annotations <- data.frame(species = list_species_hint,
                                  anno_db_hint = list_db_annotationdbi_hint,
                                  row.names = list_species_hint
                                  )


# annotation_hint() --------

#' annotation_hint ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in hint
#' @param type different interaction files provided by hint (all high-quality)
#' @param ppi_hint variable defined by ppis_hint in get_networkdata_hint()
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
#'@return ppis_hint_annotated
#'
#'@export
#'
#'
#' @examples
#' #\donttest{
#' # annotation_hint(ppi_hint, species = "HomoSapiens", version = "2024-06", type = "binary")
#' #}
#'


annotation_hint <- function(ppi_hint,
                            species,
                            version,
                            type){

# find database on corresponding species

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop("Species not found as specified by HINT,",
         "please check some valid entries by running `list_species_hint`") # stop function and print
  }

  annotation_db <-
    hint_db_annotations$anno_db_hint[match(species, hint_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_hint$Uniprot_A, ppi_hint$Uniprot_B))

    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
      entrez_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
      row.names = all_prot_ids
    )

    ppis_hint_annotated <- ppi_hint

    #adding Ensembl
    ppis_hint_annotated$Ensembl_A <-
      anno_df$ensembl_id[match(ppis_hint_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_hint_annotated$Ensembl_B <-
      anno_df$ensembl_id[match(ppis_hint_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding Entrez
    ppis_hint_annotated$Entrez_A <-
      anno_df$entrez_id[match(ppis_hint_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_hint_annotated$Entrez_B <-
      anno_df$entrez_id[match(ppis_hint_annotated$Uniprot_B, anno_df$uniprot_id)]

    return(ppis_hint_annotated)
  }

  if (is.na(annotation_db)) {
    return(ppi_hint)
  }

}


#output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B






