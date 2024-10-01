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
