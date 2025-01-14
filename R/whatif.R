# Handling MITAB files ------

# e.g. for irefindex columns (other example is mint)
# #GeneSymbol
# ppi_irefindex_df_annotated$GeneSymbol_A <- str_extract(ppi_irefindex_df_annotated$`Aliases for A`, "uniprotkb:([^\\(]+)\\(gene name\\)")
# ppi_irefindex_df_annotated$GeneSymbol_A <- gsub("uniprotkb:", "", ppi_irefindex_df_annotated$GeneSymbol_A)
# ppi_irefindex_df_annotated$GeneSymbol_A <- gsub("\\(gene name\\)", "", ppi_irefindex_df_annotated$GeneSymbol_A)
#
# ppi_irefindex_df_annotated$GeneSymbol_B <- str_extract(ppi_irefindex_df_annotated$`Aliases for B`, "uniprotkb:([^\\(]+)\\(gene name\\)")
# ppi_irefindex_df_annotated$GeneSymbol_B <- gsub("uniprotkb:", "", ppi_irefindex_df_annotated$GeneSymbol_B)
# ppi_irefindex_df_annotated$GeneSymbol_B <- gsub("\\(gene name\\)", "", ppi_irefindex_df_annotated$GeneSymbol_B)
#
# #Ensembl
# ppi_irefindex_df_annotated$Ensembl_A <- strsplit(ppi_irefindex_df_annotated$`Alternative identifier for interactor A`, "\\|") # split entry by "|" -> you get a list for each row
# ppi_irefindex_df_annotated$Ensembl_A <- lapply(ppi_irefindex_df_annotated$Ensembl_A, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
#   ensembl_entry <- grep("^ensembl:", x, value = TRUE) #search for the entry (x) in every list that starts (^) with "ensembl"
#   ensembl_entry <- gsub("ensembl:", "", ensembl_entry) #remove prefix "ensembl" by using gsub()
#   return(ensembl_entry)
# })
#
# ppi_irefindex_df_annotated$Ensembl_B <- strsplit(ppi_irefindex_df_annotated$`Alternative identifier for interactor B`, "\\|") # split entry by "|" -> you get a list for each row
# ppi_irefindex_df_annotated$Ensembl_B <- lapply(ppi_irefindex_df_annotated$Ensembl_B, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
#   ensembl_entry <- grep("^ensembl:", x, value = TRUE) #search for the entry (x) in every list that starts (^) with "ensembl"
#   ensembl_entry <- gsub("ensembl:", "", ensembl_entry) #remove prefix "ensembl" by using gsub()
#   return(ensembl_entry)
# })


# Handling MITAB files annotation (Reactome) --------------

# if (add_annotation) {
#
#
#   url_reactome_annotation_info <- "https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt"
#
#   rname <- paste0(
#     "reactome_",
#     version,
#     "_",
#     species,
#     "_annotation_info"
#   )
#
#   if (cache) {
#     # tries to fetch from the cache
#     message("Trying to fetch from cache...")
#     network_file_annotation <- fetch_NetworkHub(rname)
#   }
#
#   if (!cache | is.null(network_file_annotation)) {
#     # retrieves the file for the first time
#     message("Downloading to cache...")
#     network_file_annotation <- cache_NetworkHub(
#       rname = rname,
#       fpath = url_reactome_annotation_info
#     )
#   }
#
#   reactome_annotation_info <- vroom::vroom(network_file_annotation)
#
#   #rename columns Uniprot
#   colnames(reactome_annotation_info)[colnames(reactome_annotation_info) == "# Interactor 1 uniprot id"] <- "Uniprot_A"
#   colnames(reactome_annotation_info)[colnames(reactome_annotation_info) == "Interactor 2 uniprot id"] <- "Uniprot_B"
#
#   #merge dataframes into new dataframe
#   ppis_reactome_filtered_annotated <- merge(ppis_reactome_filtered, reactome_annotation_info, by = c("Uniprot_A", "Uniprot_B"), all.x = TRUE)
#   ppis_reactome_filtered_annotated <- ppis_reactome_filtered_annotated[
#     !(ppis_reactome_filtered_annotated$`Taxid interactor A` == "-" &
#         ppis_reactome_filtered_annotated$`Taxid interactor B` == "-"),
#   ]
#
#   #rename columns Entrez
#   colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 1 Entrez Gene id"] <- "Entrez_A"
#   colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 2 Entrez Gene id"] <- "Entrez_B"
#
#   #rename columns Ensembl
#   colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 1 Ensembl gene id"] <- "Ensembl_A"
#   colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 2 Ensembl gene id"] <- "Ensembl_B"
#
#   #extract UniProt_id
#   ppis_reactome_filtered_annotated$Uniprot_A <- str_extract(ppis_reactome_filtered_annotated$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
#   ppis_reactome_filtered_annotated$Uniprot_A <- gsub("uniprotkb:", "", ppis_reactome_filtered_annotated$Uniprot_A)
#
#   ppis_reactome_filtered_annotated$Uniprot_B <- str_extract(ppis_reactome_filtered_annotated$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
#   ppis_reactome_filtered_annotated$Uniprot_B <- gsub("uniprotkb:", "", ppis_reactome_filtered_annotated$Uniprot_B)
#
#
#   return(ppis_reactome_filtered_annotated)
# }
#
# if (!add_annotation) {
#   return(ppis_reactome_filtered)
# }
# }



# species NCBI listed in ncbi_taxid_host_organism --------------

#just a list that connects names with numbers , we only get mouse and human ppi data from innatedb

# #info_species_innatedb <- list("Mus musculus" = "10090",
#                               "Homo sapiens" = "9606",
#                               "Saccharomyces cerevisiae" ="4932",
#                               "taxid:0",
#                               "Canis lupus familiaris" = "9615",
#                               "Spodoptera frugiperda" = "7108",
#                               "Rattus norvegicus" = "10116",
#                               "Chlorocebus aethiops" = "9534",
#                               "Cricetulus griseus" = "10029",
#                               "Gallus gallus" = "9031",
#                               "Platyrrhini" = "9479",
#                               "Oryctolagus cuniculus" = "9986",
#                               "Cricetinae" = "10026",
#                               "Coturnix japonica" = "93934",
#                               "Bos taurus" = "9913",
#                               "Escherichia phage T7" = "10760",
#                               "Yeast two-hybrid vector pC-ACT.2" = "111296",
#                               "unidentified baculovirus" = "10469",
#                               "Chlorocebus aethiops aethiops" = "101841",
#                               "Escherichia coli"  ="562",
#                               "Neogale vison" = "452646",
#                               "Drosophila melanogaster" = "7227",
#                               "Hylobates lar" = "9580",
#                               "Spodoptera frugiperda multiple nucleopolyhedrovirus" = "10455",
#                               "Mesocricetus auratus" = "10036",
#                               "taxid:-1",
#                               "Cercopithecidae" = "9527",
#                               "Mustela lutreola" =  "9666",
#                               "Lepidoptera" = "7088")

# CPDB -----------------
#get_networkdata_cpdb()

#' #' get_networkdata_cpdb()
#' #'
#' #' @param species  from which species does the data come from (default human because currently only human data provided from cpdb)
#' #' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' #' @param ... 	further arguments passed to or from other methods
#' #'
#' #' @return ppis_cpdb
#' #'
#' #' @importFrom vroom vroom
#' #' @export
#' #'
#' #' @examples
#' #' \dontrun{
#' #' db_cpdb_df <- get_networkdata_cpdb(
#' #'   species = "human"
#' #' )
#' #'
#' #' db_cpdb_df
#' #' }
#'
#' get_networkdata_cpdb <- function( species = "human",
#'                                   cache = TRUE,
#'                                   ...) {
#'
#'
#'
#'   # list species is actualized for version cpdb "2021-05"
#'   # UPDATEVERSION
#'
#'   # check that the value for species is listed in cpdb
#'
#'   if (!(species %in% list_species_cpdb)) { # if species is not in the list
#'     stop("Species not found as specified by cpdb,",
#'          "please check some valid entries by running `list_species_cpdb`") # stop function and print
#'   }
#'
#'   # buildup of the resource location for the version and all
#'   ## elegantly done in another smaller utility function
#'
#'   rname <- paste0(
#'     "cpdb_",
#'     species
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
#'     # buildup from "base" cpdb url
#'     cpdb_url <-
#'       urlmaker_cpdb(
#'         species = species
#'       ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_cpdb
#'
#'     # and cache_NetworkHub to cache the file from the url source
#'     network_file <- cache_NetworkHub(
#'       rname = rname,
#'       fpath = cpdb_url
#'     )
#'   }
#'
#'   # read in the resource, whether cached or freshly downloaded
#'   ppis_cpdb <- vroom::vroom(network_file, delim = "\t", skip = 1)
#'   #ppis_cpdb <- head(read.delim(network_file, sep = " "))
#'   #
#'   message(dim(ppis_cpdb))
#'
#'   #Rename (Annotation)
#'   colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__uniprot_id"] <- "Uniprot"
#'   colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__genename"] <- "GeneSymbol"
#'   colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__entrez_gene"] <- "Entrez"
#'   colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__ensembl_gene"] <- "Ensembl"
#'
#'   return(ppis_cpdb)
#'
#' }
#'
#' ## outside of function
#'
#' list_species_cpdb <- c("human",
#'                        "mouse",
#'                        "yeast")
