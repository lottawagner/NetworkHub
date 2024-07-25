# Here we collect all the functions in order to work with the data. ----

## First function is ... -----


## What Fede has in here ...
get_networkdata_STRINGDB <- function(species,
                                     version,
                                     cache = TRUE,
                                     remap_identifiers = TRUE,
                                     remap_to = c(
                                       "ENSEMBL", "gene_symbol"
                                       # , "ENTREZ"
                                     )){

}



get_accessoryinfo_STRINGDB <- function(species,
                                       version,
                                       cache = TRUE) {

}



create_annotation_from_STRINGDBaccessory <- function(accessory_info) {
  accessory_info_df <- data.frame(
    protein_id = unique(accessory_info$`#string_protein_id`),
    row.names = unique(accessory_info$`#string_protein_id`)
  ){

  }



build_graph_STRINGDB <- function(graph_data,
                                   output_format = c("igraph", "graphnel", "sif"),
                                   min_score_threshold = NULL,
                                   ## alternatives? native for cytoscape?
                                   subset_nodes = NULL) {

  }
