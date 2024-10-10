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


#Handling MITAB files annotation (Reactome) --------------

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
