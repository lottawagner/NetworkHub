# get_networkdata_irefindex() -----------

#' get_networkdata_irefindex()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in irefindex
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with four columns (Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppi_irefindex
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_irefindex_df <- get_networkdata_irefindex(
#'   species = "Homo sapiens",
#'   version = "08-28-2023")
#'
#' db_irefindex_df
#'

get_networkdata_irefindex <- function(species,
                                 version = "08-28-2023",
                                 cache = TRUE,
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
        species = species_id,
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

  # match the annotation db with the corresponding species
  annotation_db <-
    irefindex_db_annotations$anno_db_irefindex[match(species, irefindex_db_annotations$species_irefindex)]

  if (add_annotation) {
    ppi_irefindex_df_annotated <- annotation_irefindex(ppi_irefindex = ppis_irefindex,
                                           species = species,
                                           version = version)

    return(ppi_irefindex_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_irefindex)
  #}
  }
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

species_id_irefindex <- c(
  "9606",
  "10090",
  "559292",
  "562",
  "10116",
  "4932",
  "7227",
  "6239"
)


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

# annotation_irefindex() --------

#' annotation_irefindex ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in irefindex
#' @param ppi_irefindex variable defined by ppis_irefindex in get_networkdata_irefindex()
#'
#' @importFrom AnnotationDbi mapIds
#' @import "org.Hs.eg.db"
#' @import "org.Mm.eg.db"
#' @import "org.Sc.sgd.db"
#' @import "org.EcK12.eg.db"
#' @import "org.Rn.eg.db"
#' @import "org.Sc.sgd.db"
#' @import "org.Dm.eg.db"
#' @import "org.Ce.eg.db"
#'
#'
#' @return ppis_irefindex_annotated
#'
#' @export
#'
#'
#' @examples
#'
#' annotation_irefindex(ppi_irefindex, species = "Homo sapiens", version = "08-28-2023")



annotation_irefindex <- function(ppi_irefindex,
                                 species,
                                 version) {


  # find database for corresponding species

  annotation_db <- irefindex_db_annotations$anno_db_irefindex[match(species, irefindex_db_annotations$species_irefindex)]

  ppis_irefindex_annotated <- ppi_irefindex
  all_prot_ids <- unique(c(ppi_irefindex$Uniprot_A, ppi_irefindex$Uniprot_B))

  anno_df <- data.frame(
    uniprot_id = all_prot_ids,
    ensembl_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
    entrez_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
    row.names = all_prot_ids
  )

  #adding Ensembl
  ppis_irefindex_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_irefindex_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_irefindex_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_irefindex_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppis_irefindex_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_irefindex_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_irefindex_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_irefindex_annotated$Uniprot_B, anno_df$uniprot_id)]

  return(ppis_irefindex_annotated)


}


