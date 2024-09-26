#' # get_networkdata_mint() -----------
#'
#' #' get_networkdata_mint()
#' #'
#' #' @param species  from which species does the data come from
#' #' @param version version of the data files in mint
#' #' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' #' @param ... 	further arguments passed to or from other methods
#' #'
#' #' @return ppis_mint
#' #'
#' #' @importFrom vroom vroom
#' #' @export
#' #'
#' #' @examples
#' #'
#' #' db_mint_df <- get_networkdata_mint(
#' #'   species = "human",
#' #'   version = "2021-05"
#' #' )
#' #'
#' #' db_mint_df
#' #'
#'
#' get_networkdata_mint <- function( species,
#'                                  version = "2021-05",
#'                                  cache = TRUE,
#'                                  add_annotation = TRUE,
#'                                  ...) {
#'
#'
#'
#'   # list species is actualized for version mint "2021-05"
#'   # UPDATEVERSION
#'
#'   # check that the value for species is listed in mint
#'
#'   if (!(species %in% list_species_mint)) { # if species is not in the list
#'     stop("Species not found as specified by mint,",
#'          "please check some valid entries by running `list_species_mint`") # stop function and print
#'   }
#'
#'   # buildup of the resource location for the version and all
#'   ## elegantly done in another smaller utility function
#'
#'   rname <- paste0(
#'     "mint_",
#'     species,
#'     "_v",
#'     version
#'   ) # definition of the resource name
#'
#'   if (cache) {
#'     # tries to fetch from the cache
#'     message("Trying to fetch from cache...")
#'     network_file <- fetch_NetworkHub(rname)
#'   } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file
#'
#'   if (!cache | is.null(network_file)) {
#'     # retrieves the file for the first time
#'     message("Downloading to cache...")
#'     # buildup from "base" mint url
#'     mint_url <-
#'       urlmaker_mint(
#'         species = species,
#'         version = version
#'       ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_mint
#'
#'     # and cache_NetworkHub to cache the file from the url source
#'     network_file <- cache_NetworkHub(
#'       rname = rname,
#'       fpath = mint_url
#'     )
#'   }
#'
#'   # read in the resource, whether cached or freshly downloaded
#'   ppis_mint <- vroom::vroom(network_file)
#'   #ppis_mint <- head(read.delim(network_file, sep = " "))
#'   #
#'   message(dim(ppis_mint))
#'
#'   #Uniprot
#'   colnames(ppis_mint)[colnames(ppis_mint) == "uniprot1"] <- "Uniprot_A"
#'   colnames(ppis_mint)[colnames(ppis_mint) == "uniprot2"] <- "Uniprot_B"
#'
#'   #GeneSymbol
#'   colnames(ppis_mint)[colnames(ppis_mint) == "symbol1"] <- "GeneSymbol_A"
#'   colnames(ppis_mint)[colnames(ppis_mint) == "symbol2"] <- "GeneSymbol_B"
#'
#'   annotation_db <-
#'     mint_db_annotations$anno_db_mint[match(species, mint_db_annotations$species)]
#'
#'   if (add_annotation && !is.na(annotation_db)) {
#'     ppi_mint_df_annotated <- annotation_mint(ppi_mint = ppis_mint,
#'                                            species = species,
#'                                            version = version)
#'     return(ppi_mint_df_annotated)
#'   }
#'
#'   if (add_annotation && is.na(annotation_db)) {
#'     message("Annotation database for the species is not implemented yet.\n",
#'             "Next time define add_annotation in get_networkdata_mint(..., add_annotation = FALSE, ...)\n",
#'             "You will get ppis_mint only containing annotation columns for Uniprot_A/B & GeneSymbol_A/B.")
#'     return(ppis_mint)
#'   }
#'
#'   if (!add_annotation) {
#'     return(ppis_mint)
#'   }
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # outside of function ----------
#'
#' list_species_mint <- c ("alpaca",
#'                        "cat",
#'                        "chicken",
#'                        "cow",
#'                        "dog",
#'                        "duck",
#'                        "fly",
#'                        "guinea_pig",
#'                        "horse",
#'                        "human",
#'                        "mouse",
#'                        "pig",
#'                        "rabbit",
#'                        "rat",
#'                        "sheep",
#'                        "turkey",
#'                        "worm",
#'                        "yeast")
#'
#' list_db_annotationdbi_mint <- c(NA,
#'                                NA,
#'                                "org.Gg.eg.db",
#'                                "org.Bt.eg.db",
#'                                "org.Cf.eg.db",
#'                                NA,
#'                                "org.Dm.eg.db",
#'                                NA,
#'                                NA,
#'                                "org.Hs.eg.db",
#'                                "org.Mm.eg.db",
#'                                "org.Ss.eg.db",
#'                                NA,
#'                                "org.Rn.eg.db",
#'                                NA,
#'                                NA,
#'                                "org.Ce.eg.db",
#'                                "org.Sc.sgd.db"
#' )
#'
#'
#' mint_db_annotations <- data.frame(species = list_species_mint,
#'                                  anno_db_mint = list_db_annotationdbi_mint,
#'                                  row.names = list_species_mint
#' )
#'
#'
#' # annotation_mint() --------
#'
#' #' annotation_mint ()
#' #'
#' #' @param species  from which species does the data come from
#' #' @param version version of the data files in mint
#' #' @param type different interaction files provided by mint (all high-quality)
#' #' @param ppi_mint variable defined by ppis_mint in get_networkdata_mint()
#' #'
#' #' @importFrom AnnotationDbi mapIds
#' #' @import org.Gg.eg.db
#' #' @import org.Bt.eg.db
#' #' @import org.Cf.eg.db
#' #' @import org.Dm.eg.db
#' #' @import org.Hs.eg.db
#' #' @import org.Mm.eg.db
#' #' @import org.Ss.eg.db
#' #' @import org.Rn.eg.db
#' #' @import org.Ce.eg.db
#' #' @import org.Sc.sgd.db
#' #'
#' #'
#' #'@return ppis_mint_annotated
#' #'
#' #'@export
#' #'
#' #'
#' #' @examples
#' #'
#' #' # annotation_mint(ppi_mint, species = "HomoSapiens", version = "2024-06", type = "binary")
#' #' #TODO: what can I do here as ppi_funcoup is not defined in annotation_mint()?
#'
#'
#' annotation_mint <- function(ppi_mint,
#'                            species,
#'                            version) {
#'   # find database on corresponding species
#'
#'   if (!(species %in% list_species_mint)) { # if species is not in the list
#'     stop("Species not found as specified by mint,",
#'          "please check some valid entries by running `list_species_mint`") # stop function and print
#'   }
#'
#'   annotation_db <-
#'     mint_db_annotations$anno_db_mint[match(species, mint_db_annotations$species)]
#'
#'   if (!is.na(annotation_db)) {
#'     all_prot_ids <- unique(c(ppi_mint$Uniprot_A, ppi_mint$Uniprot_B))
#'     anno_df <- data.frame(
#'       uniprot_id = all_prot_ids,
#'       ensembl_id = mapIds(
#'         get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
#'       entrez_id = mapIds(
#'         get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
#'       row.names = all_prot_ids
#'     )
#'
#'
#'     ppis_mint_annotated <- ppi_mint
#'
#'     #adding Ensembl
#'     ppis_mint_annotated$Ensembl_A <-
#'       anno_df$ensembl_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]
#'     ppis_mint_annotated$Ensembl_B <-
#'       anno_df$ensembl_id[match(ppis_mint_annotated$Uniprot_B, anno_df$uniprot_id)]
#'
#'     #adding Entrez
#'     ppis_mint_annotated$Entrez_A <-
#'       anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_A, anno_df$uniprot_id)]
#'     ppis_mint_annotated$Entrez_B <-
#'       anno_df$entrez_id[match(ppis_mint_annotated$Uniprot_B, anno_df$uniprot_id)]
#'
#'     return(ppis_mint_annotated)
#'   }
#'
#'
#'   if (is.na(annotation_db)) {
#'     return(ppi_mint)
#'   }
#' }
#'
#' #output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B
#'
#'
#'
#'
#'
#'
#'
