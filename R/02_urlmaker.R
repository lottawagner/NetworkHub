# STRINGDB -----

## stringdb species -----

#' info_species_stringdb
#'
#' @param version #current version of the data files in stringdb
#'
#' @return
#' @export
#'
#' @examples
#'

info_species_stringdb <- function(version = "12.0"){

  # use sprintf function to insert current version into the url
  url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt", version)

  # read.delim the data from the species text file (columns separated using a delimiter)
  df_species <- read.delim(url(url_species_stringdb))

  # return the datafile to the caller
  return(df_species)
}


## stringdb urlmaker -----

#' urlmaker_stringdb
#'
#' @param type # type of information that is stored in the file
#' @param species # from which species does the data come from
#' @param version
#'
#' @return
#' @export
#'
#' @examples
#'

urlmaker_stringdb <- function(type = "PPI", # or protein info for stringdb
                              species = "Homo sapiens",
                              version = "12.0") {

  # "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"

  #match the information about the species with the corresponding data file (PPI/protein info)

  info_species <- info_species_stringdb(version = version)
  species_id <- info_species$X.taxon_id[
    match(species, info_species$official_name_NCBI)
  ]

  if (type == "PPI") {
    url <- sprintf(
      "https://stringdb-downloads.org/download/protein.links.v%s/%s.protein.links.v%s.txt.gz",
      version,
      species_id,
      version
    )
  } else if (type == "protein_info") {
    url <- sprintf(
      "https://stringdb-downloads.org/download/protein.info.v%s/%s.protein.info.v%s.txt.gz",
      version,
      species_id,
      version
    )
  }

  return(url)
}


