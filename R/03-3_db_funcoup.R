# get_networkdata_funcoup() -----------

#' get_networkdata_funcoup()
#'
#' @param species  from which species does the data come from c( "A.thaliana", "B.subtilis", "B.taurus", "C.elegans","C.familiaris", "C.intestinalis", "D.melanogatser", "D.rerio", "E.coli", "G.gallus", "H.sapiens", "M.jannaschii", "M.musculus", "O.sativa", "P.falciparum", "R.norvegicus", "S.cerevisae", "S.pombe", "S.scrofa", "S.solfataricus")
#' @param version version of the data files in funcoup
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with four columns ( Entrez_ID and Ensembl_ID )
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_funcoup
#' @return ppis_annotated_funcoup
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_funcoup_df <- get_networkdata_funcoup(
#'   species = "H.sapiens",
#'   version = "5.0",
#' )
#'
#' db_funcoup_df
#'

get_networkdata_funcoup <- function(species = "H.sapiens",
                                 version = "5.0",
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {



  # list species is actualized for version funcoup v.5.0 build 2020-09
  # UPDATEVERSION

  # check that the value for species is listed in funcoup

  if (!(species %in% list_species_funcoup)) { # if species is not in the list
    stop("Species not found as specified by FunCoup,",
         "please check some valid entries by running `list_species_funcoup`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "funcoup_",
    species,
    "_v",
    version
  ) # definition of the resource name

  #TODO problem with caching what is happening? it is always downloading again (3x)
  if (cache) {
    # tries to fetch from the cache
    message("Fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("File not in cache, downloading to cache...")
    # buildup from funcoup url
    funcoup_url <-
      urlmaker_funcoup(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_funcoup

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = funcoup_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_funcoup <- vroom::vroom(network_file)
  #ppis_funcoup <- head(read.delim(network_file, sep = "\t"))

  message(dim(ppis_funcoup))
  #colnames(ppis_funcoup)

  colnames(ppis_funcoup)[colnames(ppis_funcoup) == "2:Gene1"] <- "Ensembl_A"
  colnames(ppis_funcoup)[colnames(ppis_funcoup) == "3:Gene2"] <- "Ensembl_B"

  head(ppis_funcoup)

  if (add_annotation) {
    ppis_funcoup_df_annotated <- annotation_funcoup(ppis_funcoup = ppis_funcoup,
                                             species = species,
                                             version = version)
    message(dim(ppis_funcoup_df_annotated))
    return(ppis_funcoup_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_funcoup)
    message(dim(ppis_funcoup))
  }


}


# outside of function ----------

#TODO do I need to load the library?
BiocManager::install("AnnotationDbi", force = TRUE)
# BiocManager::install("org.At.tair.db")
# BiocManager::install("org.Bt.eg.db")
# BiocManager::install("org.Ce.eg.db")
# BiocManager::install("org.Cf.eg.db")
# BiocManager::install("org.Dm.eg.db")
# BiocManager::install("org.Dr.eg.db")
# BiocManager::install("org.EcK12.eg.db")
# BiocManager::install("org.Gg.eg.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
# BiocManager::install("org.Pf.plasmo.db")
# BiocManager::install("org.Rn.eg.db")
# BiocManager::install("org.Sc.sgd.db")
# BiocManager::install("org.Ss.eg.db")


library(AnnotationDbi)
# library(org.At.tair.db)
# library(org.Bt.eg.db)
# library(org.Ce.eg.db)
# library(org.Cf.eg.db)
# library(org.Dm.eg.db)
# library(org.Dr.eg.db)
# library(org.EcK12.eg.db)
# library(org.Gg.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
# library(org.Pf.plasmo.db)
# library(org.Rn.eg.db)
# library(org.Sc.sgd.db)
# library(org.Ss.eg.db)


list_species_funcoup <- c( "A.thaliana",
                           "B.subtilis",
                           "B.taurus",
                           "C.elegans",
                           "C.familiaris",
                           "C.intestinalis",
                           "D.melanogatser",
                           "D.rerio",
                           "E.coli",
                           "G.gallus",
                           "H.sapiens",
                           "M.jannaschii",
                           "M.musculus",
                           "O.sativa",
                           "P.falciparum",
                           "R.norvegicus",
                           "S.cerevisae",
                           "S.pombe",
                           "S.scrofa",
                           "S.solfataricus")

list_db_annotationdbi_funcoup <- c("org.At.tair.db",
                                   NA,
                                   "org.Bt.eg.db",
                                   "org.Ce.eg.db",
                                   "org.Cf.eg.db",
                                   NA,
                                   "org.Dm.eg.db",
                                   "org.Dr.eg.db",
                                   "org.EcK12.eg.db",
                                   "org.Gg.eg.db",
                                   "org.Hs.eg.db",
                                   NA,
                                   "org.Mm.eg.db",
                                   NA,
                                   "org.Pf.plasmo.db",
                                   "org.Rn.eg.db",
                                   "org.Sc.sgd.db",
                                   NA,
                                   "org.Ss.eg.db",
                                   NA)


funcoup_db_annotations <- data.frame(
  species = list_species_funcoup,
  anno_db_funcoup = list_db_annotationdbi_funcoup,
  row.names = list_species_funcoup
)


# annotation_funcoup() --------

#' annotation_funcoup ()
#'
#' @param species from which species does the data come from c( "A.thaliana", "B.subtilis", "B.taurus", "C.elegans","C.familiaris", "C.intestinalis", "D.melanogatser", "D.rerio", "E.coli", "G.gallus", "H.sapiens", "M.jannaschii", "M.musculus", "O.sativa", "P.falciparum", "R.norvegicus", "S.cerevisae", "S.pombe", "S.scrofa", "S.solfataricus")
#' @param version version of the data files in funcoup
#' @param ppis_funcoup variable defined by ppis_funcoup in get_networkdata_funcoup()
#'
#'@return ppis_funcoup_annotated
#' @export
#'
#' @examples
#'
#' annotation_funcoup(ppis_funcoup, species = "H.sapiens", version = "5.0")
#'


annotation_funcoup <- function(ppis_funcoup,
                            species,
                            version){

  # find database on corresponding species

  if (!(species %in% list_species_funcoup)) { # if species is not in the list
    stop("Species not found as specified by FunCoup,",
         "please check some valid entries by running `list_species_funcoup`") # stop function and print
  }

  if (species %in% list_species_funcoup) {
    annotation_db <-
      funcoup_db_annotations$anno_db_funcoup[match(species, funcoup_db_annotations$species)]

    if (is.na(annotation_db)) {
      stop("Annotation database for the species is not implemented yet.")
    }
  }

  all_gene_ids <- unique(c(ppis_funcoup$Ensembl_A, ppis_funcoup$Ensembl_B))

  print(dim(all_gene_ids))

  anno_df <- data.frame(
    ensembl_id = all_gene_ids,
    gene_symbol = mapIds(
      get(annotation_db), keys = all_gene_ids, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"),
    uniprot_id = mapIds(
      get(annotation_db), keys = all_gene_ids, keytype = "ENSEMBL", column = "UNIPROT", multiVals = "first"),
    entrez_id = mapIds(
      get(annotation_db), keys = all_gene_ids, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first"),
    row.names = all_gene_ids
  )

  ppis_funcoup_annotated <- ppis_funcoup

  ppis_funcoup_annotated$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppis_funcoup_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_funcoup_annotated$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppis_funcoup_annotated$Ensembl_B, anno_df$ensembl_id)]

  ppis_funcoup_annotated$Uniprot_A <-
    anno_df$uniprot_id[match(ppis_funcoup_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_funcoup_annotated$Uniprot_B <-
    anno_df$uniprot_id[match(ppis_funcoup_annotated$Ensembl_B, anno_df$ensembl_id)]

  ppis_funcoup_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_funcoup_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_funcoup_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_funcoup_annotated$Ensembl_B, anno_df$ensembl_id)]

  #TODO maybe create a dataframe that only contains 8 columns (GeneSymbol_A/B, uniprot A/B, Ensmebl A/B, Entrez A/B)?
  return(ppis_funcoup_annotated)

}






