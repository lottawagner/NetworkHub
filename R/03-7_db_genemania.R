# get_networkdata_genemania() -----------

#' get_networkdata_genemania()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in genemania
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with four columns (Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_genemania
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_genemania_df <- get_networkdata_genemania(
#'   species = "Homo_sapiens",
#'   version = "current"
#' )
#'
#' db_genemania_df
#' }

get_networkdata_genemania <- function( species = "Homo_sapiens",
                                       version = "current",
                                       cache = TRUE,
                                       add_annotation = TRUE,
                                       ...) {



  # list species is actualized for version genemania "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in genemania

  if (!(species %in% list_species_genemania)) { # if species is not in the list
    stop("Species not found as specified by genemania,",
         "please check some valid entries by running `list_species_genemania`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "genemania_",
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
    # buildup from "base" genemania url
    genemania_url <-
      urlmaker_genemania(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_genemania

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = genemania_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_genemania <- vroom::vroom(network_file)
  #ppis_genemania <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_genemania))

  #Uniprot
  colnames(ppis_genemania)[colnames(ppis_genemania) == "Gene_A"] <- "Uniprot_A"
  colnames(ppis_genemania)[colnames(ppis_genemania) == "Gene_B"] <- "Uniprot_B"

  annotation_db <-
    genemania_db_annotations$anno_db_genemania[match(species, genemania_db_annotations$species)]

  if (add_annotation) {
    ppi_genemania_df_annotated <- annotation_genemania(ppi_genemania = ppis_genemania,
                                           species = species,
                                           version = version)
    return(ppi_genemania_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_genemania)
  }
}

# outside of function ----------

list_species_genemania <- c("Arabidopsis_thaliana",
                            "Caenorhabditis_elegans",
                            "Danio_rerio",
                            "Drosophila_melanogaster",
                            "Escherichia_coli",
                            "Homo_sapiens",
                            "Mus_musculus",
                            "Rattus_norvegicus",
                            "Saccharomyces_cerevisiae")

list_db_annotationdbi_genemania <- c("org.At.tair.db",
                                     "org.Ce.eg.db",
                                     "org.Cf.eg.db",
                                     "org.Dm.eg.db",
                                     "org.EcK12.eg.db",
                                     "org.Hs.eg.db",
                                     "org.Mm.eg.db",
                                     "org.Rn.eg.db",
                                     "org.Sc.sgd.db")


genemania_db_annotations <- data.frame(species = list_species_genemania,
                                 anno_db_genemania = list_db_annotationdbi_genemania,
                                 row.names = list_species_genemania
)


# annotation_genemania() --------

#' annotation_genemania ()
#'
#' @param ppi_genemania variable defined by ppis_genemania in get_networkdata_genemania()
#' @param species  from which species does the data come from
#' @param version version of the data files in genemania
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.At.tair.db
#' @import org.Ce.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.EcK12.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#'
#' @return ppis_genemania_annotated
#'
#' @export
#'
#' @examples
#' #\dontrun{
#' # annotation_genemania <- annotation_genemania(ppi_genemania, species = "Homo_sapiens", version = "current")
#' #annotation_genemania
#' #}



annotation_genemania <- function(ppi_genemania,
                           species,
                           version) {
  # find database on corresponding species

  if (!(species %in% list_species_genemania)) { # if species is not in the list
    stop("Species not found as specified by genemania,",
         "please check some valid entries by running `list_species_genemania`") # stop function and print
  }

  annotation_db <-
    genemania_db_annotations$anno_db_genemania[match(species, genemania_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_genemania$Uniprot_A, ppi_genemania$Uniprot_B))
    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      gene_symbol = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL"),
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL"),
      entrez_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID"),
      row.names = all_prot_ids
    )


    ppis_genemania_annotated <- ppi_genemania

    #adding GeneSymbol
    ppis_genemania_annotated$GeneSymbol_A <-
      anno_df$gene_symbol[match(ppis_genemania_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_genemania_annotated$GeneSymbol_B <-
      anno_df$gene_symbol[match(ppis_genemania_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding Ensembl
    ppis_genemania_annotated$Ensembl_A <-
      anno_df$ensembl_id[match(ppis_genemania_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_genemania_annotated$Ensembl_B <-
      anno_df$ensembl_id[match(ppis_genemania_annotated$Uniprot_B, anno_df$uniprot_id)]

    #adding Entrez
    ppis_genemania_annotated$Entrez_A <-
      anno_df$entrez_id[match(ppis_genemania_annotated$Uniprot_A, anno_df$uniprot_id)]
    ppis_genemania_annotated$Entrez_B <-
      anno_df$entrez_id[match(ppis_genemania_annotated$Uniprot_B, anno_df$uniprot_id)]

    return(ppis_genemania_annotated)
  }


  if (is.na(annotation_db)) {
    return(ppi_genemania)
  }
}

#output: dataframe containing 4 columns:  Uniprot_A  Uniprot_B Gene_A Gene_B







