

# get_networkdata_reactome() -----------

#' get_networkdata_reactome()
#'
#' @param species  default value = "taxid:9606(Homo sapiens)" - from which species does the data come from
#' @param version version of the data files in reactome
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with six columns (GeneSymbol,Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_reactome
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' db_reactome_df <- get_networkdata_reactome(species = "taxid:9606(Homo sapiens)",
#'                                            version = "current"
#'                                            )
#' db_reactome_df
#'

get_networkdata_reactome <- function(species = "taxid:9606(Homo sapiens)",
                                     version = "current",
                                     cache = TRUE,
                                     add_annotation = TRUE,
                                     ...) {

  # list species is actualized for version reactome "current" (sept 2024)
  # UPDATEVERSION

  # check that the value for species is listed in reactome

  if (!(species %in% list_species_reactome)) { # if species is not in the list
    stop("Species not found as specified by reactome,",
         "please check some valid entries by running `list_species_reactome`") # stop function and print
  }


  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "reactome_v_",
    version
    )
  #UPDATEVERSION

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" reactome url
    reactome_url <-
      urlmaker_reactome(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_reactome

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = reactome_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_reactome <- vroom::vroom(network_file)
  #ppis_reactome <- head(read.delim(network_file, sep = " "))

  # If you want to check all species contained in the file
  #species_A <- unique(ppis_reactome$`Taxid interactor A`)
  #species_B <- unique(ppis_reactome$`Taxid interactor B`)
  #info_species_reactome <- union(species_A, species_B)


  #Filter for species

  ppis_reactome_filtered <- ppis_reactome[(ppis_reactome$`Taxid interactor A` == species ) &
                                            (ppis_reactome$`Taxid interactor B` == species ),
  ]

  # rename columns

  #Uniprot
  colnames(ppis_reactome_filtered)[colnames(ppis_reactome_filtered) == "#ID(s) interactor A"] <- "Uniprot_A"
  colnames(ppis_reactome_filtered)[colnames(ppis_reactome_filtered) == "ID(s) interactor B"] <- "Uniprot_B"

  #extract Uniprot_id
  ppis_reactome_filtered$Uniprot_A <- str_extract(ppis_reactome_filtered$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
  ppis_reactome_filtered$Uniprot_A <- gsub("uniprotkb:", "", ppis_reactome_filtered$Uniprot_A)

  ppis_reactome_filtered$Uniprot_B <- str_extract(ppis_reactome_filtered$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
  ppis_reactome_filtered$Uniprot_B <- gsub("uniprotkb:", "", ppis_reactome_filtered$Uniprot_B)



  if (add_annotation) {

    ppi_reactome_filtered_annotated <- annotation_reactome(ppi_reactome = ppis_reactome_filtered,
                                                           species = species,
                                                           version = version)
    return(ppi_reactome_filtered_annotated)
  }

  if (!add_annotation) {
    return(ppis_reactome_filtered)
  }
}


# Outside function ----------


#SPECIES NetworkHub uses from reactome
list_species_reactome <- c("taxid:9913(Bos taurus)",
                           "taxid:6239(Caenorhabditis elegans)",
                           "taxid:9615(Canis familiaris)",
                           "taxid:7227(Drosophila melanogaster)",
                           "taxid:562(Escherichia coli)",
                           "taxid:9031(Gallus gallus)",
                           "taxid:9606(Homo sapiens)",
                           "taxid:10090(Mus musculus)",
                           "taxid:10116(Rattus norvegicus",
                           "taxid:4932(Saccharomyces cerevisiae)",
                           "taxid:8355(Xenopus laevis)"
                          )


list_db_annotationdbi_reactome <- c("org.Bt.eg.db",
                                   "org.Ce.eg.db",
                                   "org.Cf.eg.db",
                                   "org.Dm.eg.db",
                                   "org.EcK12.eg.db",
                                   "org.Gg.eg.db",
                                   "org.Hs.eg.db",
                                   "org.Mm.eg.db",
                                   "org.Rn.eg.db",
                                   "org.Sc.sgd.db",
                                   "org.Xl.eg.db"
                                   )


reactome_db_annotations <- data.frame(species = list_species_reactome,
                                 anno_db_reactome = list_db_annotationdbi_reactome,
                                 row.names = list_species_reactome
)


# annotation_reactome() --------

#' annotation_reactome ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in reactome
#' @param ppi_reactome variable defined by ppis_reactome in get_networkdata_reactome()
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Bt.eg.db
#' @import org.Ce.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.EcK12.eg.db
#' @import org.Gg.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#' @import org.Xl.eg.db
#'
#'@return ppis_reactome_annotated
#'
#'@export
#'
#'
#' @examples
#' # annotation_reactome(ppi_reactome, species = "taxid:9606(Homo sapiens)", version = "current")

annotation_reactome <- function(ppi_reactome,
                                species,
                                version)
                                {
  # find database on corresponding species

  if (!(species %in% list_species_reactome)) { # if species is not in the list
    stop("Species not found as specified by reactome,",
         "please check some valid entries by running `list_species_reactome`") # stop function and print
  }

  annotation_db <-
    reactome_db_annotations$anno_db_reactome[match(species, reactome_db_annotations$species)]


  all_prot_ids <- unique(na.omit(c(ppi_reactome$Uniprot_A, ppi_reactome$Uniprot_B)))

  anno_df <- data.frame(uniprot_id = all_prot_ids,
                        gene_symbol <- mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "UNIPROT",
                                              column = "SYMBOL"
                                              ),
                        ensembl_id <- mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "UNIPROT",
                                              column = "ENSEMBL"
                                            ),
                        entrez_ids <- mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "UNIPROT",
                                             column = "ENTREZID"
                                            ),
                        row.names = all_prot_ids
                        )

  ppis_reactome_annotated <- ppi_reactome

  #adding GeneSymbol
  ppis_reactome_annotated$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppis_reactome_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_reactome_annotated$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppis_reactome_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Ensembl
  ppis_reactome_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_reactome_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_reactome_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_reactome_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppis_reactome_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_reactome_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_reactome_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_reactome_annotated$Uniprot_B, anno_df$uniprot_id)]

  return(ppis_reactome_annotated)
}










