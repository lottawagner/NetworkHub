# get_networkdata_mint() -----------

#' get_networkdata_mint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in mint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_mint
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_mint_df <- get_networkdata_mint(
#'   species = "Homo Sapiens",
#'   version = "current"
#' )
#'
#' db_mint_df
#'

get_networkdata_mint <- function(species,
                                 version = "current",
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {



  # list species is actualized for version mint "current"
  # UPDATEVERSION

  # check that the value for species is listed in mint

  if (!(species %in% list_species_mint)) { # if species is not in the list
    stop("Species not found as specified by mint,",
         "please check some valid entries by running `list_species_mint`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "mint_",
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
    # buildup from "base" mint url
    mint_url <-
      urlmaker_mint(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_mint

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = mint_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_mint <- vroom::vroom(network_file)
  #ppis_mint <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_mint))

  colnames(ppis_mint) <- c("Identifier_A", "Identifier_B", "Alternative identifier for interactor A", "Alternative identifier for interactor B", "Aliases for A", "Aliases for B", "method", "First author", "pubmed", "NCBI Taxonomy identifier for interactor A","NCBI Taxonomy identifier for interactor B", "interaction_type", "Source databases and identifiers", "Interaction identifier(s) in the corresponding source database", "Confidence score")
  # ppis_mint$Uniprot_A <- sub("uniprotkb:", "", ppis_mint$Identifier_A)
  # ppis_mint$Uniprot_A [grepl("intact:", ppis_mint$Identifier_A)] <- NA
  # ppis_mint$Uniprot_A [grepl("ensembl", ppis_mint$Identifier_A)] <- NA
  # ppis_mint$Uniprot_A [grepl("chebi:", ppis_mint$Identifier_A)] <- NA
  # ppis_mint$Uniprot_B <- sub("uniprotkb:", "", ppis_mint$Identifier_B)
  # ppis_mint$Uniprot_B [grepl("intact:", ppis_mint$Identifier_B)] <- NA
  # ppis_mint$Uniprot_B [grepl("ensembl", ppis_mint$Identifier_B)] <- NA
  # ppis_mint$Uniprot_B [grepl("chebi:", ppis_mint$Identifier_B)] <- NA
  #
  # annotation_db <-
  #   mint_db_annotations$anno_db_mint[match(species, mint_db_annotations$species)]


  ppis_mint$GeneSymbol_A <- sub("...:", "","(gene..." ppis_mint$`Aliases for A`)
  ppis_mint_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]

  if (add_annotation) {
    ppi_mint_df_annotated <- annotation_mint(ppi_mint = ppis_mint,
                                           species = species,
                                           version = version)
    return(ppi_mint_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_mint)
  }
}


# outside of function ----------

list_species_mint <- c ("Homo Sapiens",
                        "Mus Musculus",
                        "Drosophila Melanogaster",
                        "Saccharomyces Cerevisiae")

list_db_annotationdbi_mint <- c("org.Hs.eg.db",
                                "org.Mm.eg.db",
                                "org.Dm.eg.db",
                                "org.Sc.sgd.db"
                                )


mint_db_annotations <- data.frame(species = list_species_mint,
                                 anno_db_mint = list_db_annotationdbi_mint,
                                 row.names = list_species_mint
)


# annotation_mint() --------

#' annotation_mint ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in mint
#' @param type different interaction files provided by mint (all high-quality)
#' @param ppi_mint variable defined by ppis_mint in get_networkdata_mint()
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Sc.sgd.db
#'
#'
#'@return ppis_mint_annotated
#'
#'@export
#'
#'
#' @examples
#'
#' # annotation_mint(ppi_mint, species = "HomoSapiens", version = "2024-06", type = "binary")
#' #TODO: what can I do here as ppi_funcoup is not defined in annotation_mint()?


# extract_annotation <- function(gene_symbol,
#                                ensembl,
#                                entrez,
#                                uniprot) {
#
#   if (gene_symbol == ){
#     parts <- strsplit(gene_symbol, "\\|")[[1]]
#     for (part in parts) {
#       if (grepl("uniprotkb:", part) && grepl("\\(gene name\\)", part)) {
#         return(sub(".*uniprotkb:([^\\(]+).*", "\\1", part))
#       }
#     }
#     return(NA)
# }

annotation_mint <- function(ppi_mint,
                           species,
                           version) {
  # find database on corresponding species

  if (!(species %in% list_species_mint)) { # if species is not in the list
    stop("Species not found as specified by mint,",
         "please check some valid entries by running `list_species_mint`") # stop function and print
  }

  annotation_db <-
    mint_db_annotations$anno_db_mint[match(species, mint_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_mint$Uniprot_A, ppi_mint$Uniprot_B))
    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
      entrez_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
      gene_symbol = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL"),
      row.names = all_prot_ids
    )


    ppis_mint_annotated <- ppi_mint

    #adding Ensembl
    ppis_mint_annotated$Ensembl_A <-
      anno_df$ensembl_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_mint_annotated$Ensembl_B <-
      anno_df$ensembl_id[match(ppis_mint_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding Entrez
    ppis_mint_annotated$Entrez_A <-
      anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_mint_annotated$Entrez_B <-
      anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding GeneSymbol
    ppis_mint_annotated$GeneSymbol_A <-
      anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_mint_annotated$GeneSymbol_B <-
      anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_B, anno_df$uniprot_id)]

    return(ppis_mint_annotated)
  }


  if (is.na(annotation_db)) {
    return(ppi_mint)
  }
}

#output:







