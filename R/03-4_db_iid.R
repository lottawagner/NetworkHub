

# get_networkdata_iid() -----------

#' get_networkdata_iid()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in iid
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_iid
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_iid_df <- get_networkdata_iid(
#'   species = "human",
#'   version = "2021-05"
#' )
#'
#' db_iid_df
#'

get_networkdata_iid <- function( species,
                                 version = "2021-05",
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {



  # list species is actualized for version iid "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in iid

  if (!(species %in% list_species_iid)) { # if species is not in the list
    stop("Species not found as specified by iid,",
         "please check some valid entries by running `list_species_iid`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "iid_",
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
    # buildup from "base" iid url
    iid_url <-
      urlmaker_iid(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_iid

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = iid_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_iid <- vroom::vroom(network_file)
  #ppis_iid <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_iid))

  #Uniprot
  colnames(ppis_iid)[colnames(ppis_iid) == "uniprot1"] <- "Uniprot_A"
  colnames(ppis_iid)[colnames(ppis_iid) == "uniprot2"] <- "Uniprot_B"

  #GeneSymbol
  colnames(ppis_iid)[colnames(ppis_iid) == "symbol1"] <- "GeneSymbol_A"
  colnames(ppis_iid)[colnames(ppis_iid) == "symbol2"] <- "GeneSymbol_B"

  annotation_db <-
      iid_db_annotations$anno_db_iid[match(species, iid_db_annotations$species)]

  if (add_annotation && !is.na(annotation_db)) {
    ppi_iid_df_annotated <- annotation_iid(ppi_iid = ppis_iid,
                                             species = species,
                                             version = version)
    return(ppi_iid_df_annotated)
  }

  if (add_annotation && is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_iid(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_iid only containing annotation columns for Uniprot_A/B & GeneSymbol_A/B.")
    return(ppis_iid)
  }

  if (!add_annotation) {
      return(ppis_iid)
  }
}

# outside of function ----------

list_species_iid <- c ("alpaca",
                       "cat",
                       "chicken",
                       "cow",
                       "dog",
                       "duck",
                       "fly",
                       "guinea_pig",
                       "horse",
                       "human",
                       "mouse",
                       "pig",
                       "rabbit",
                       "rat",
                       "sheep",
                       "turkey",
                       "worm",
                       "yeast")

list_db_annotationdbi_iid <- c(NA,
                               NA,
                               "org.Gg.eg.db",
                               "org.Bt.eg.db",
                               "org.Cf.eg.db",
                               NA,
                               "org.Dm.eg.db",
                               NA,
                               NA,
                               "org.Hs.eg.db",
                               "org.Mm.eg.db",
                               "org.Ss.eg.db",
                               NA,
                               "org.Rn.eg.db",
                               NA,
                               NA,
                               "org.Ce.eg.db",
                               "org.Sc.sgd.db"
                              )


iid_db_annotations <- data.frame(species = list_species_iid,
                                  anno_db_iid = list_db_annotationdbi_iid,
                                  row.names = list_species_iid
)


# annotation_iid() --------

#' annotation_iid ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in iid
#' @param type different interaction files provided by iid (all high-quality)
#' @param ppi_iid variable defined by ppis_iid in get_networkdata_iid()
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Gg.eg.db
#' @import org.Bt.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Ss.eg.db
#' @import org.Rn.eg.db
#' @import org.Ce.eg.db
#' @import org.Sc.sgd.db
#'
#'
#'@return ppis_iid_annotated
#'
#'@export
#'
#'
#' @examples
#'
#' # annotation_iid(ppi_iid, species = "HomoSapiens", version = "2024-06", type = "binary")
#' #TODO: what can I do here as ppi_funcoup is not defined in annotation_iid()?


annotation_iid <- function(ppi_iid,
                            species,
                            version) {
  # find database on corresponding species

  if (!(species %in% list_species_iid)) { # if species is not in the list
    stop("Species not found as specified by iid,",
         "please check some valid entries by running `list_species_iid`") # stop function and print
  }

  annotation_db <-
      iid_db_annotations$anno_db_iid[match(species, iid_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_iid$Uniprot_A, ppi_iid$Uniprot_B))
    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
      entrez_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
      row.names = all_prot_ids
    )


    ppis_iid_annotated <- ppi_iid

    #adding Ensembl
    ppis_iid_annotated$Ensembl_A <-
      anno_df$ensembl_id[match(ppis_iid_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_iid_annotated$Ensembl_B <-
      anno_df$ensembl_id[match(ppis_iid_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding Entrez
    ppis_iid_annotated$Entrez_A <-
      anno_df$entrez_id[match(ppis_iid_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_iid_annotated$Entrez_B <-
      anno_df$entrez_id[match(ppis_iid_annotated$Uniprot_B, anno_df$uniprot_id)]

    return(ppis_iid_annotated)
  }


  if (is.na(annotation_db)) {
    return(ppi_iid)
  }
}

#output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B







