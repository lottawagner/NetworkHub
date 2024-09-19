#' get_networkdata_hint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in hint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param ... 	further arguments passed to or from other methods
#' @param type different interaction files provided by hint (all high-quality)
#'
#' @return ppis_hint
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_hint_df <- get_networkdata_hint(
#'   species = "Homo sapiens",
#'   version = "2024-06",
#'   type = "binary"
#' )
#'
#' db_hint_df
#'

get_networkdata_hint <- function(species,
                                 version,
                                 type = "binary",
                                 cache = TRUE,
                                 ...) {

  # matching for the species...
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

  # list species is actualized for version HINT "2024-06"
  # UPDATEVERSION

  # check that the value for species is listed in HINT

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop("Species not found as specified by HINT,",
         "please check some valid entries by running `list_species_hint`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))

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

  #add annotation

  if (species == "HomoSapiens") {
    library(org.Hs.eg.db)
    annotation_db <- org.Hs.eg.db
  } else {
    stop("Annotation database for the species is not implemented yet.")
  }

  protein_ids_A <- unique(ppis_hint$Uniprot_A)
  protein_ids_B <- unique(ppis_hint$Uniprot_B)

  # Annotations for Uniprot_A (EntrezID, EnsemblID)
  annotations_A <- AnnotationDbi::select(
    annotation_db,
    keys = protein_ids_A,
    columns = c("ENTREZID", "ENSEMBL"),
    keytype = "UNIPROT"
  )

  # Annotations for Uniprot_A (EntrezID, EnsemblID)
  annotations_B <- AnnotationDbi::select(
    annotation_db,
    keys = protein_ids_B,
    columns = c("ENTREZID", "ENSEMBL"),
    keytype = "UNIPROT"
  )

  # define the names vor the columns
  colnames(annotations_A) <- c("Uniprot_A", "Entrezid_A", "Ensembl_A")
  colnames(annotations_B) <- c("Uniprot_B", "Entrezid_B", "Ensembl_B")

  # merge with ppis_hint
  ppis_hint_annotated <- merge(ppis_hint, annotations_A, by = "Uniprot_A", all.x = TRUE)
  ppis_hint_annotated <- merge(ppis_hint_annotated, annotations_B, by = "Uniprot_B", all.x = TRUE)

  message(dim(ppis_hint_annotated)[1])

  # return ppis-hint
  return(ppis_hint)
}

## annotation

annotation_hint <- function(annotation_info){

  annotation_info_df <- data.frame(
    protein_id = unique(annotation_info$`Gene_A`),
    row.names = unique(annotation_info$`Gene_B`)
  )

}


#output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B



