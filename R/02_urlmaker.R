
# What's inside ?
# 1. For each database one/more functions that help to identify the variables in the corresponding url_maker function
# 2. For each database one function to create an url to download/cache the data from the db
# -> Why? because using an url where you have to define the version, species and some other variables, the user can easily access the data without the need of google search

# STRINGDB -----

# For stringdb we can use a file that defines the names of the organisms and correlating species_id to tell the url_maker function what to put inside the url by chosing a name





#' urlmaker_stringdb()
#'
#' @param type type of information that is stored in the file
#' @param species from which species does the data come from
#' @param version version of the data files in stringdb
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_stringdb(type = "PPI",
#'                     species = "Homo sapiens",
#'                     version = "12.0")
urlmaker_stringdb <- function(type = "PPI",
                              species = "Homo sapiens",
                              version = "12.0") {

  type <- match.arg(type, c("PPI", "protein_info")) # define the possible types in StringDB
  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  info_species <- info_species_stringdb(version = version) # assign info_species to the info about stringdb species of the chosen version using function info_species_stringdb() defined in 03-01_db_stringdb.R

  # check that the value for species is listed in HINT
  if (!(species %in% info_species$official_name_NCBI)) {
    stop("Species not found as specified by STRINGDB, ",
         "please check some valid entries by running `info_species_stringdb()`")
  }

  # match the information about the species (id, name) with the corresponding data file (PPI/protein info)

  species_id <- info_species$X.taxon_id[
    match(species, info_species$official_name_NCBI)
  ]

  # here we define the url for the different type options

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

  # return value is "url"
  return(url)
}

# HINT ------------------------------------------

#links:
# binary: https://hint.yulab.org/download-raw/2024-06/HomoSapiens_binary_hq.txt
# cocomplex: https://hint.yulab.org/download-raw/2024-06/HomoSapiens_cocomp_hq.txt
# https://hint.yulab.org/download-raw/2024-06/SaccharomycesCerevisiae_binary_hq.txt



#' urlmaker_hint()
#'
#' @param type  interaction types in HINT, default value = "binary"

# "binary" = binary
# "cocomp" = co-complex
# "lcb" = literature curated binary
# "lcc" = literature curated co-complex

#' @param species types listed in list_species_hint depending on current version, default value = "HomoSapiens" #UPDATEVERSION
#' @param version version of HINT, always written as year - month ( 2024-06, 2020-08, ...)
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_hint(type = "binary",
#'               species = "HomoSapiens",
#'               version = "2024-06")
urlmaker_hint <- function(type = "binary", # default value for type = "binary"
                          species = "HomoSapiens", # default value for species #UPDATEVERSION
                          version) { # version is always written as year - month ( 2024-06, 2020-08, ...) #UPDATEVERSION

  # avoid errors
  stopifnot(is.character(species))  # make sure to type in a species name as character
  stopifnot(is.character(version))  # make sure to type in a version as character
  stopifnot(length(version) == 1)   # make sure to type in a version with the length == 1


  # define the opportunities for species in HINT

  list_species_hint <- c("HomoSapiens",
                         "SaccharomycesCerevisiae",
                         "SchizosaccharomycesPombe",
                         "MusMusculus",
                         "DrosophilaMelanogaster",
                         "CaenorhabditisElegans",
                         "ArabidopsisThaliana",
                         "EscherichiaColi",
                         "RattusNorvegicus",
                         "OryzaSativa")

                          # list species is actualized for version HINT "2024-06"
                          # UPDATEVERSION

  # check that the value for species is listed in HINT

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop("Species not found as specified by HINT,",
         "please check some valid entries by running `list_species_hint`") # stop function and print
  }

  # define the opportunities for the four interaction types in HINT:

  # "binary" = binary
  # "cocomp" = co-complex
  # "lcb" = literature curated binary
  # "lcc" = literature curated co-complex

  type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))



  # create the url depending on the type, version and species

  url <- sprintf("https://hint.yulab.org/download-raw/%s/%s_%s_hq.txt",
                 version,
                 species,
                 type)


  # return value is "url"
  return(url)

}





# FunCoup ----------


#' urlmaker_funcoup()
#'
#' @param version version of FunCoup , default value = "5.0", value as type "5.0" #UPDATEVERSION
#' @param species types listed in list_species_funcoup depending on current version, default value = "H.sapiens" #UPDATEVERSION
#'
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_funcoup( version = "4.1",
#'                   species = "B.taurus")
urlmaker_funcoup <- function(version = "5.0", # default value = "5.0", value as type "5.0"
                             species = "H.sapiens") { # default value = "H.sapiens", value as type first letter (capital) of first name, full second name with small letters

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1)  # make sure to type in a version with the length == 1

  # define the opportunities for species in FunCoup

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
                          "S.solfataricus"
                          )

                          # list species is actualized for version FC5.0
                          # UPDATEVERSION

  # check that the value for species is listed in FunCoup
  if (!(species %in% list_species_funcoup)) {
    stop("Species not found as specified by FunCoup,",
         "please check some valid entries and version of `list_species_funcoup`")
  }

  # using sprintf to define the url of FunCoup by looking at the corresponding version and species

  # current version is stored on website in download folder #UPDATEVERSION
  if (version == "5.0") {
    url <- sprintf("https://funcoup.org/downloads/download.action?type=network&instanceID=24480085&fileName=FC%s_%s_full.gz",
                 version,
                 species)

  # archived versions are in archive folder #UPDATEVERSION
  } else if (version != "5.0") {
    url <- sprintf("https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_full.gz",
                   version,
                   version,
                   species)
  }

  # return value is "url
  return(url)
}






# ONLY HUMAN
# IID --------------------

#' urlmaker_iid()
#'
#' @param species types listed in list_species_iid depending on current version, default value = "human" #UPDATEVERSION
#' @param version version of IID , default value = "2021-05", value as type "2021-05" #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_iid(species = "mouse")
urlmaker_iid <- function(species = "human", #
                         version = "2021-05") { # version of IID not updated since

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1 ) # make sure to type in a version with the length == 1


  # define the opportunities for species in IID

  list_species_iid <- c ("alpaca",
                         "cat",
                         "chicken",
                         "cow",
                         "dog",
                         "duck",
                         "fly",
                         "guinea_pig",
                         "horse",
                         "human",
                         "mouse",
                         "pig",
                         "rabbit",
                         "rat",
                         "sheep",
                         "turkey",
                         "worm",
                         "yeast")
                          # list species is actualized for version IID version 2021-05
                          # UPDATEVERSION

  # check that the value for species is listed in IID
  if (!species %in% list_species_iid){
    stop("Species not found as specified by IID,",
         "please check some valid entries of `list_species_iid` and on the webiste 'https://iid.ophid.utoronto.ca/search_by_proteins/'")
  }

  # create the url depending on the species
  url <- sprintf("https://iid.ophid.utoronto.ca/static/download/%s_annotated_PPIs.txt.gz",
  species)

  #return the url for IID and the corresponding species
  return(url)
}


# CPDB - ONLY HUMAN -----------------------------


#' urlmaker_cpdb()
#'
#' @param species currently only human (default value)
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_cpdb(species = "human")
#'
urlmaker_cpdb <- function (species = "human") { #default value = human because at the moment (05.08.2024) only human

  list_species_cpdb <- c("human", "mouse", "yeast")

  stopifnot(is.character(species)) # make sure to type in a species name as character

  url <- sprintf("http://cpdb.molgen.mpg.de/download/ConsensusPathDB_%s_PPI.gz",
                 species)
  return(url)

}
# HuRi - ONLY HUMAN ----------------

#' urlmaker_huri
#'
#' @param species default value = "human", because this database only provides human data
#' @param type differnet datasets , more information on "http://www.interactome-atlas.org/about/"
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#'
#' urlmaker_huri(species = "human",
#'              type = "HI-union")
urlmaker_huri <- function (species = "human", # default value human, because this database only provides human data
                           type = "HI-union") {

  list_species_huri <- c("human")

  stopifnot(is.character(species)) # make sure to type in a species name as character

  if (!(species %in% list_species_huri)) { # if species is not in the list
    stop("Species not found as specified by HuRi,",
         "please check some valid entries by running `list_species_huri`") # stop function and print
  }

  # the datafile "HI-union" is an aggregate of all PPIs identified in
  # HI-I-05,
  # HI-II-14
  # HuRI
  # Venkatesan-09
  # Yu-11
  # Yang-16
  # Test space screens-19

  type <- c("HI-union")

  url <- sprintf("http://www.interactome-atlas.org/data/%s.tsv",
                 type)
  return(url)

}



# TODO


# TODO
# BIOGRID - TODO --------------------------------------


#https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.235/BIOGRID-ORGANISM-4.4.235.tab3.zip
#https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.235/BIOGRID-ORGANISM-4.4.235.tab3.zip


#' info_species_biogrid
#'
#' @param version version of the data files in stringdb
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' # TODO
# info_species_zip_biogrid <- function

#' urlmaker_biogrid

urlmaker_biogrid <- function(type = "PPI", # or protein info for stringdb
                              species = "Homo sapiens",
                              version = "4.4.235") {

}







# IntAct - TODO -----

species_intact <- c("human",
                    "Escherichia",
                    "mouse"
                    )

#
# urlmaker_intact <- function(type = "PPI", # could be also protein info # QUESTION: WHY?
#                             species = "human",
#                             version = "current")
#   https://ftp.ebi.ac.uk/pub/databases/intact/%s/%s/species

# https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/

#Human Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/human.zip

#E.Coli Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/Escherichia.zip

#Mouse Zip -> https://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/species/mouse.zip


# CORUM - TODO ------------------------------
#Idee: all datei herunterladen und in get_networkdata_corum definieren welchen organimus man anschauen will
# https://mips.helmholtz-muenchen.de/corum/download/releases/current/humanComplexes.txt.zip
# https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip
# GeneMania - TODO-------------

# https://pages.genemania.org/data/ TODO: There are folders containing
# multiple files for each species ... how can I work on that?
# combine the single file manually?





