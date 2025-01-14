# get_networkdata_cpdb() -----------

#' get_networkdata_cpdb()
#'
#' @param species  from which species does the data come from (default human because currently only human data provided from cpdb)
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_cpdb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_cpdb_df <- get_networkdata_cpdb(
#'   species = "human"
#' )
#'
#' db_cpdb_df
#' }

get_networkdata_cpdb <- function(species = "human",
                                 cache = TRUE,
                                 ...) {



  # list species is actualized for version cpdb "2021-05"
  # UPDATEVERSION

  # check that the value for species is listed in cpdb

  if (!(species %in% list_species_cpdb)) { # if species is not in the list
    stop("Species not found as specified by cpdb,",
         "please check some valid entries by running `list_species_cpdb`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "cpdb_",
    species
  ) # definition of the resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" cpdb url
    cpdb_url <-
      urlmaker_cpdb(
        species = species
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_cpdb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = cpdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_cpdb <- vroom::vroom(network_file, delim = "\t", skip = 1)
  #ppis_cpdb <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_cpdb))

  #Rename (Annotation)
  colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__uniprot_id"] <- "Uniprot"
  colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__genename"] <- "GeneSymbol"
  colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__entrez_gene"] <- "Entrez"
  colnames(ppis_cpdb)[colnames(ppis_cpdb) == "interaction_participants__ensembl_gene"] <- "Ensembl"

  return(ppis_cpdb)

}

# outside of function ----------

list_species_cpdb <- c("human",
                       "mouse",
                       "yeast")
