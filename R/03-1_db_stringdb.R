# Here we collect all the functions in order to work with the data. ----

##
get_networkdata_STRINGDB <- function(species,
                                     version,
                                     cache = TRUE,
                                     remap_identifiers = TRUE,
                                     remap_to = c(
                                       "ENSEMBL", "gene_symbol"
                                       # , "ENTREZ"
                                     ),
                                     ...){

  # matching for the species...
  species_id <- NULL # looking for the corresponding id in the following

  info_species <- info_species_stringdb(version = version) # we find the information about species_name and species_id in the file info_species_stringdb
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
    match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]


  # if (is.null(species_id))
  #   stop("No match found!")

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function



  rname <- paste0(
    "STRINGDB_",
    species,
    "_v",
    version
  )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- get_NetworkHub(rname)
  }

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Retrieving to cache...")
    # buildup from "base" STRINGDB url
    stringdb_url <-
      urlmaker_stringdb(
        type = "PPI",
        species = species,
        version = version
      )


    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = stringdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  sdb_ppis <- vroom::vroom(network_file)
  # sdb_ppis <- read.delim(network_file, sep = " ")

  # WE COULD STOP HERE, but we won't - at the end...


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
