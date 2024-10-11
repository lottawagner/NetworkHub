# get_networkdata_pathwaycommons() -----------

#' get_networkdata_pathwaycommons()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in pathwaycommons
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with four columns (Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_pathwaycommons
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \donttest{
#' db_pathwaycommons_df <- get_networkdata_pathwaycommons(
#'   species = "human",
#'   version = "v12"
#' )
#'
#' db_pathwaycommons_df
#' }

get_networkdata_pathwaycommons <- function( species = "human",
                                       version = "v12",
                                       cache = TRUE,
                                       add_annotation = TRUE,
                                       ...) {



  # list species is actualized for version pathwaycommons "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in pathwaycommons

  if (!(species %in% list_species_pathwaycommons)) { # if species is not in the list
    stop("Species not found as specified by pathwaycommons,",
         "please check some valid entries by running `list_species_pathwaycommons`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "pathwaycommons_",
    species,
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
    # buildup from "base" pathwaycommons url
    pathwaycommons_url <-
      urlmaker_pathwaycommons(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_pathwaycommons

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = pathwaycommons_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_pathwaycommons <- vroom::vroom(network_file)
  #ppis_pathwaycommons <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_pathwaycommons))

  #GeneSymbol
  colnames(ppis_pathwaycommons)[colnames(ppis_pathwaycommons) == "PARTICIPANT_A"] <- "GeneSymbol_A"
  colnames(ppis_pathwaycommons)[colnames(ppis_pathwaycommons) == "PARTICIPANT_B"] <- "GeneSymbol_B"

  annotation_db <-
    pathwaycommons_db_annotations$anno_db_pathwaycommons[match(species, pathwaycommons_db_annotations$species)]

  if (add_annotation) {
    ppi_pathwaycommons_df_annotated <- annotation_pathwaycommons(ppi_pathwaycommons = ppis_pathwaycommons,
                                                       species = species,
                                                       version = version)
    return(ppi_pathwaycommons_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_pathwaycommons)
  }
}

# outside of function ----------

list_species_pathwaycommons <- c("human")

list_db_annotationdbi_pathwaycommons <- c("org.Hs.eg.db")


pathwaycommons_db_annotations <- data.frame(species = list_species_pathwaycommons,
                                       anno_db_pathwaycommons = list_db_annotationdbi_pathwaycommons,
                                       row.names = list_species_pathwaycommons
)


# annotation_pathwaycommons() --------

#' annotation_pathwaycommons ()
#'
#' @param ppi_pathwaycommons variable defined by ppis_pathwaycommons in get_networkdata_pathwaycommons()
#' @param species from which species does the data come from
#' @param version version of the data files in pathwaycommons
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db

#'
#' @return ppis_pathwaycommons_annotated
#'
#' @export
#'
#' @examples
#' #\donttest{
#' #annotation_pathwaycommons <- annotation_pathwaycommons(ppi_pathwaycommons, species = "human", version = "v12")
#' #annotation_pathwaycommons
#' #}



annotation_pathwaycommons <- function(ppi_pathwaycommons,
                                      species,
                                      version) {
  # find database on corresponding species

  if (!(species %in% list_species_pathwaycommons)) { # if species is not in the list
    stop("Species not found as specified by pathwaycommons,",
         "please check some valid entries by running `list_species_pathwaycommons`") # stop function and print
  }

  annotation_db <-
    pathwaycommons_db_annotations$anno_db_pathwaycommons[match(species, pathwaycommons_db_annotations$species)]

  all_prot_ids <- unique(c(ppi_pathwaycommons$GeneSymbol_A, ppi_pathwaycommons$GeneSymbol_B))
  anno_df <- data.frame(
    genesymbol = all_prot_ids,
    uniprot_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "UNIPROT"),
    ensembl_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "ENSEMBL"),
    entrez_id = mapIds(
      get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "ENTREZID"),
    row.names = all_prot_ids
    )


  ppis_pathwaycommons_annotated <- ppi_pathwaycommons

  #adding Uniprot
  ppis_pathwaycommons_annotated$Uniprot_A <-
    anno_df$uniprot_id[match(ppis_pathwaycommons_annotated$GeneSymbol_A, anno_df$genesymbol)]
  ppis_pathwaycommons_annotated$Uniprot_B <-
    anno_df$uniprot_id[match(ppis_pathwaycommons_annotated$GeneSymbol_B, anno_df$genesymbol)]

  #adding Ensembl
  ppis_pathwaycommons_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_pathwaycommons_annotated$GeneSymbol_A, anno_df$genesymbol)]
  ppis_pathwaycommons_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_pathwaycommons_annotated$GeneSymbol_B, anno_df$genesymbol)]

  #adding Entrez
  ppis_pathwaycommons_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_pathwaycommons_annotated$GeneSymbol_A, anno_df$genesymbol)]
  ppis_pathwaycommons_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_pathwaycommons_annotated$GeneSymbol_B, anno_df$genesymbol)]

  return(ppis_pathwaycommons_annotated)

}









