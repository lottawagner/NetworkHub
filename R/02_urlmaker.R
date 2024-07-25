
# What's inside ?
# 1. For each database one/more functions that help to identify the variables in the corresponding url_maker function
# 2. For each database one function to create an url to download/cache the data from the db
# -> Why? because using an url where you have to define the version, species and some other variables, the user can easily access the data without the need of google search


# STRINGDB -----

# For stringdb we can use a file that defines the names of the organisms and correlating species_id to tell the url_maker function what to put inside the url by chosing a name

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

  #match the information about the species (id, name) with the corresponding data file (PPI/protein info)

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

# BIOGRID ----

https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.235/BIOGRID-ORGANISM-4.4.235.tab3.zip
https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.235/BIOGRID-ORGANISM-4.4.235.tab3.zip

## biogrid species -----

#' info_species_biogrid
#'
#' @param version #current version of the data files in stringdb
#'
#' @return
#' @export
#'
#' @examples
#'

info_species_zip_biogrid <- function
## biogrid urlmaker -----

#' urlmaker_biogrid

urlmaker_biogrid <- function(type = "PPI", # or protein info for stringdb
                              species = "Homo sapiens",
                              version = "4.4.235") {

}







# IntAct -----

Human Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/human.zip

E.Coli Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/Escherichia.zip

Mouse Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/mouse.zip



