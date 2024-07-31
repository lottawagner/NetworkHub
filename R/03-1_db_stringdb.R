# Here we collect all the functions in order to work with the data. ----

## Functions we use in the get_networkdata_stringdb ----

get_accessoryinfo_stringdb <- function(species,
                                       version,
                                       cache = TRUE) {

  # species is matched to the id in the info file

  species_id <- NULL

  info_species <- info_species_stringdb(version = version) # we find the information about species_name and species_id in the file info_species_stringdb
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
    match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]

  # rname of accessory info defined

  rname <- paste0(
    "STRINGDB_accessoryInfo_",
    species,
    "v_",
    version
  )

  # fetch (& cache) the accessory info_file

  ## fetch from the cache
  if (cache) {
    message("Trying to fetch from cache...")
    proteininfo_file <- fetch_NetworkHub(rname)
  }

  ## cache file and fetch from the cache

  if (!cache | is.null(proteininfo_file)) #QUESTIONFEDE: Why do we need to say both?
   message("Retrieving to cache...")
   stringdb_url <- urlmaker_stringdb(
     type = "protein_info",
     species = species,
     version = version
   )
   proteininfo_file <- cache_NetworkHub(
     rname = rname,
     fpath = stringdb_url
   )
 }

# reading in the fetched information file
proteininfo_stringdb <- vroom::vroom(proteininfo_file)

return(proteininfo_stringdb)
}

## Annotation ----

#QUESTIONFEDE: Is the idea here to be able to tell the user what gene/protein we are looking at in different annotation types?
#If so, where does the accessoryinfo comes from?

create_annotation_from_stringdbaccessory <- function(accessory_info) {
  accessory_info_df <- data.frame(
    protein_id = unique(accessory_info$`#string_protein_id`),
    row.names = unique(acessory_info$`#string_protein_id`)
  )

  # use NA to assign a previous "value" to the variables
  accessory_info_df$ensemble_id <- NA
  accessory_info_df$gene_symbol <- NA
  # accessory_info_df$entrez_id <- NA

  df_ensembl <- accessory_info[accessory_info$source == "Ensembl_gene",]
}

## Graph ----

build_graph_STRINGDB <- function(graph_data,
                                   output_format = c("igraph", "graphnel", "sif"),
                                   min_score_threshold = NULL,
                                   ## alternatives? native for cytoscape?
                                   subset_nodes = NULL) {

  }


## final function to get the data from stringdb in the way we want it ----


get_networkdata_stringdb <- function(species,
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


  if (is.null(species_id))
    stop("No match found!")

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "STRINGDB_",
    species,
    "_v",
    version
  ) # definition of the  resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- get_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function get_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Retrieving to cache...")
    # buildup from "base" stringdb url
    stringdb_url <-
      urlmaker_stringdb(
        type = "PPI",
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_stringdb


    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = stringdb_url
    ) # and cache_NetworkHub to cache the file from the url source
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_stringdb <- vroom::vroom(network_file)
  #ppis_stringdb <- read.delim(network_file, sep = " ")

  # WE COULD STOP HERE, but we won't - at the end...

  #get_acessoryinfo_stringdb
  #annotation
  #


}




