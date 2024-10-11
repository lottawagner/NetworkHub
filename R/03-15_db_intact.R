# get_networkdata_intact() -----------

#' get_networkdata_intact()
#'
#' @param species  default value = "taxid:9606(human)|taxid:9606(Homo sapiens)" - from which species does the data come from
#' @param version default value = "current" version of the data files in intact
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with six columns (GeneSymbol,Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_intact
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \donttest{
#' db_intact_df <- get_networkdata_intact(species = "taxid:9606(human)|taxid:9606(Homo sapiens)",
#'                                            version =  "current"
#'                                            )
#' db_intact_df
#' }

get_networkdata_intact <- function(species = "taxid:9606(human)|taxid:9606(Homo sapiens)",
                                   version = "current",
                                   cache = TRUE,
                                   add_annotation = TRUE,
                                    ...) {
  # UPDATEVERSION

  # buildup of the resource location for the version

  rname <- paste0("intact_v_",
                  version
                 )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" intact url
    intact_url <-
      urlmaker_intact(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_intact

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = intact_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  message("...vroom takes some time (intact.txt ~7 GB) ...")
  ppis_intact <- vroom::vroom(unz(network_file, "intact.txt"))


  #Filter for species

  species_A <- unique(ppis_intact$`Taxid interactor A`)
  species_B <- unique(ppis_intact$`Taxid interactor B`)
  list_species_intact <- union(species_A, species_B)

  if (!(species %in% list_species_intact)) { # if species is not in the list
    stop("Species not found as specified by intact,",
         "please check some valid entries by running `list_species_intact`") # stop function and print
  }

  #filter out data according to chosen species
  ppis_intact_filtered <- ppis_intact[(ppis_intact$`Taxid interactor A` == species ) &
                                      (ppis_intact$`Taxid interactor B` == species ),]

  # rename columns
  #Uniprot
  colnames(ppis_intact_filtered)[colnames(ppis_intact_filtered) == "#ID(s) interactor A"] <- "Uniprot_A"
  colnames(ppis_intact_filtered)[colnames(ppis_intact_filtered) == "ID(s) interactor B"] <- "Uniprot_B"

  #extract Uniprot_id
  ppis_intact_filtered$Uniprot_A <- str_extract(ppis_intact_filtered$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
  ppis_intact_filtered$Uniprot_A <- gsub("uniprotkb:", "", ppis_intact_filtered$Uniprot_A)

  ppis_intact_filtered$Uniprot_B <- str_extract(ppis_intact_filtered$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
  ppis_intact_filtered$Uniprot_B <- gsub("uniprotkb:", "", ppis_intact_filtered$Uniprot_B)

  if (add_annotation) {

    if (!(species %in% list_common_species_intact)) { # if species is not in the list
      stop("Species not in `list_common_species_intact´!",
           "Annotation for this species is not provided") # stop function and print
    }

    ppis_intact_filtered_annotated <- annotation_intact(ppi_intact = ppis_intact_filtered,
                                                        species = species,
                                                        version = version)
    message("...added annotation :)")
    return(ppis_intact_filtered_annotated)
  }

  if (!add_annotation){
    return(ppis_intact_filtered)
  }

}



list_common_species_intact <- c("taxid:3702(arath)|taxid:3702(\"Arabidopsis thaliana (Mouse-ear cress)\")",
                                "taxid:9913(bovin)|taxid:9913(\"Bos taurus (Bovine)\")",
                                "taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)",
                                "taxid:9615(canlf)|taxid:9615(\"Canis familiaris (dog)\")",
                                "taxid:7227(drome)|taxid:7227(\"Drosophila melanogaster (Fruit fly)\")",
                                "taxid:83333(ecoli)|taxid:83333(\"Escherichia coli (strain K12)\")",
                                "taxid:9031(chick)|taxid:9031(\"Gallus gallus (Chicken)\")",
                                "taxid:9606(human)|taxid:9606(Homo sapiens)",
                                "taxid:10090(mouse)|taxid:10090(Mus musculus)",
                                "taxid:10116(rat)|taxid:10116(\"Rattus norvegicus (Rat)\")",
                                "taxid:559292(yeast)|taxid:559292(Saccharomyces cerevisiae)",
                                "taxid:9823(pig)|taxid:9823(\"Sus scrofa (Pig)\")",
                                "taxid:8355(xenla)|taxid:8355(\"Xenopus laevis (African clawed frog)\")"
                                )

list_db_annotationdbi_intact <- c("org.At.tair.db",
                                  "org.Bt.eg.db",
                                  "org.Ce.eg.db",
                                  "org.Cf.eg.db",
                                  "org.Dm.eg.db",
                                  "org.EcK12.eg.db",
                                  "org.Gg.eg.db",
                                  "org.Hs.eg.db",
                                  "org.Mm.eg.db",
                                  "org.Rn.eg.db",
                                  "org.Sc.sgd.db",
                                  "org.Ss.eg.db",
                                  "org.Xl.eg.db"
                                  )


intact_db_annotations <- data.frame(species = list_common_species_intact,
                                    anno_db_intact = list_db_annotationdbi_intact,
                                    row.names = list_common_species_intact
                                    )


# annotation_intact() --------

#' annotation_intact ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in intact
#' @param ppi_intact variable defined by ppis_intact in get_networkdata_intact()
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
#' @import org.At.tair.db
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
#' @import org.Ss.eg.db
#' @import org.Xl.eg.db
#'
#'@return ppis_intact_annotated
#'
#'@export
#'
#'
#' @examples
#' #\donttest{
#' #annotation_intact(ppi_intact, species = "taxid:9606(Homo sapiens)", version = "current")
#' #}

annotation_intact <- function(ppi_intact,
                              species,
                              version){

  annotation_db <-
    intact_db_annotations$anno_db_intact[match(species, intact_db_annotations$species)]

  #create a list that contains all uniprot ids (but not NA)
  all_prot_ids <- unique(na.omit(c(ppi_intact$Uniprot_A, ppi_intact$Uniprot_B)))

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

  ppis_intact_annotated <- ppi_intact

  #adding GeneSymbol
  ppis_intact_annotated$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppis_intact_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_intact_annotated$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppis_intact_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Ensembl
  ppis_intact_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_intact_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_intact_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_intact_annotated$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppis_intact_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_intact_annotated$Uniprot_A, anno_df$uniprot_id)]
  ppis_intact_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_intact_annotated$Uniprot_B, anno_df$uniprot_id)]

  return(ppis_intact_annotated)
}


