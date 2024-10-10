
# What's inside ?
# 1. For each database one/more functions that help to identify the variables in the corresponding url_maker function
# 2. For each database one function to create an url to download/cache the data from the db
# -> Why? because using an url where you have to define the version, species and some other variables, the user can easily access the data without the need of google search

# Comment on comments
#UPDATEVERSION - for databases where version is not defined in url
#SPECIESDEFINITION - for databases where species is not defined in url
#CURRENTVERSION - for databases only providing data for current version (no archive)


# STRINGDB -----

# For stringdb we can use a file that defines the names of the organisms and correlating species_id to tell the url_maker function what to put inside the url by choosing a name

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
#' url_stringdb <- urlmaker_stringdb(type = "PPI",
#'                                   species = "Homo sapiens",
#'                                   version = "12.0")
#' url_stringdb
urlmaker_stringdb <- function(type = "PPI",
                              species = "Homo sapiens",
                              version = "12.0") {

  type <- match.arg(type, c("PPI", "protein_info")) # define the possible types in StringDB
  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  info_species <- info_species_stringdb(version = version) # assign info_species to the info about stringdb species of the chosen version using function info_species_stringdb() defined in 03-01_db_stringdb.R

  # check that the value for species is listed in StringDB
  if (!(species %in% info_species$official_name_NCBI)) {
    stop("Species not found as specified by STRINGDB, ",
         "please check some valid entries by running `info_species_stringdb()`")
  }

  # match the information about the species (id, name) with the corresponding data file (PPI/protein info)

  species_id <- info_species$X.taxon_id[
    match(species, info_species$official_name_NCBI)
  ]

  # here we define the url for the different type options
  #QUESTION: my link is stringdb-downloads and yours is:
  #https://stringdb-static.org/download/protein.links.full.v%s/%s.protein.links.full.v%s.txt.gz

  if (type == "PPI") {
    url <- sprintf(
      "https://stringdb-downloads.org/download/protein.links.v%s/%s.protein.links.v%s.txt.gz",
      version,
      species_id,
      version
    )
  } else if (type == "protein_info") {
    url <- sprintf(
      "https://stringdb-downloads.org/download/protein.aliases.v%s/%s.protein.aliases.v%s.txt.gz",
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
#' url_hint <- urlmaker_hint( type = "binary",
#'                            species = "HomoSapiens",
#'                            version = "2024-06")
#'url_hint
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
#' @param type interaction types in FunCoup (for current version compact recommended, for older versions only full possible)
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_funcoup <- urlmaker_funcoup(version = "4.1",
#'                                 species = "B.taurus",
#'                                 type = "full")
#' url_funcoup
urlmaker_funcoup <- function(version = "5.0", # default value = "5.0", value as type "5.0"
                             species = "H.sapiens", # default value = "H.sapiens", value as type first letter (capital) of first name, full second name with small letters
                             type = c("compact", "full")) { #for current version compact recommended, for older versions only full possible

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1)  # make sure to type in a version with the length == 1

  # define the opportunities for species in FunCoup 5.0

  message("NOTE: FunCoup provides different species for the differnet version and uses different names. Check on <https://funcoup.org/archive/> for the correct definition.")

  # using sprintf to define the url of FunCoup by looking at the corresponding version and species

  # current version is stored on website in download folder #UPDATEVERSION

  if (version == "5.0") {

    list_species_funcoup_5.0 <- c( "A.thaliana",
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

    if (!species %in% list_species_funcoup_5.0){
      stop("Species not found as specified by Funcoup version 5.0,",
           "please check some valid entries on the webiste 'https://funcoup.org/archive/'")
    }


    url <- sprintf("https://funcoup.org/downloads/download.action?type=network&instanceID=24480085&fileName=FC%s_%s_%s.gz",
                  version,
                  species,
                  type)

  # archived versions are in archive folder #UPDATEVERSION
  }

  if ( version == "4.1" || version == "4.0"){

    list_species_funcoup_4._ <- c( "A.thaliana",
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
                                   "M.musculus",
                                   "O.sativa",
                                   "P.falciparum",
                                   "R.norvegicus",
                                   "S.cerevisae",
                                   "S.pombe"
                                  )
    if (!species %in% list_species_funcoup_4._){
      stop("Species not found as specified by Funcoup version 4.0 or 4.1,",
           "please check some valid entries on the webiste 'https://funcoup.org/archive/'")
    }

    url <- sprintf("https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_full.gz",
                   version,
                   version,
                   species)

    if ( type == "compact" ) {
      message("NOTE: For version 4.0 & 4.1: Only the full datafile is provided by funcoup.")
    }
  }


  if (version == "3.0") {

    list_species_funcoup_3.0 <- c( "A.thaliana",
                                   "C.elegans",
                                   "C.familiaris",
                                   "C.intestinalis",
                                   "D.melanogatser",
                                   "D.rerio",
                                   "G.gallus",
                                   "H.sapiens",
                                   "M.jannaschii",
                                   "M.musculus",
                                   "P.falciparum",
                                   "R.norvegicus",
                                   "S.cerevisae"
                                  )


    if (!species %in% list_species_funcoup_3.0){
        stop("Species not found as specified by Funcoup version 3.0,",
             "please check some valid entries on the webiste 'https://funcoup.org/archive/'")
    }

    url <- sprintf("https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_%s.gz",
                   version,
                   version,
                   species,
                   type)
  }


  if ( version == "2.0") {

    list_species_funcoup_2.0 <- c("athaliana",
                                  "celegans",
                                  "cfamiliaris",
                                  "cintestinalis",
                                  "dmelanogaster",
                                  "drerio",
                                  "ggallus",
                                  "hsapiens",
                                  "mmusculus",
                                  "rnorvegicus",
                                  "scerevisiae")


    if (!species %in% list_species_funcoup_2.0){
      stop("Species not found as specified by Funcoup version 2.0,",
           "please check some valid entries on the webiste 'https://funcoup.org/archive/'")
    }

    url <- sprintf("https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=%s.%s.pfc01.tsv.gz",
                   version,
                   species,
                   type)
  }

  if ( version == "1.0") {
    stop("No url provided by NetworkHub for version 1.0")
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
#' url_iid <- urlmaker_iid(species = "mouse")
#' url_iid
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



# iRefIndex ----------------------------


#' #' urlmaker_irefindex()
#'
#' @param species from which species does the data come from, default value = "Homo sapiens"
#' @param version version of data files in iRefIndex, default value = "2023-08-29" #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_irefindex <- urlmaker_irefindex(species = "Homo sapiens",
#'                                     version = "08-28-2023")
#' url_irefindex
urlmaker_irefindex <- function(species,
                               version = "08-28-2023"){ #default value for version = "08-28-2023" #UPDATEVERSION

  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  # create a list that contains species names and corresponding IDs

  list_species_irefindex <- c ( "Homo sapiens",
                                "Mus musculus",
                                "Saccharomyces cerevisiae S288C",
                                "Escherichia",
                                "Rattus norvegicus",
                                "Saccharomyces cerevisiae",
                                "Drosophila melanogaster",
                                "Caenorhabditis elegans"
  )

  # check that the value for species is listed in iRefIndex
  if (!species %in% list_species_irefindex) {
    stop("Species not found as specified by iRefIndex, ",
         "please check some valid entries in info_species_irefindex_id or on the website 'https://irefindex.vib.be/wiki/index.php/README_MITAB2.6_for_iRefIndex_20.0#Column_number:_29_.28Host_organism_taxid.29`")
  }

  # fetch the species_id from info_species_irefindex
  species_row <- irefindex_db_annotations[irefindex_db_annotations$species_irefindex == species, ]
  species_id <- species_row$species_id

  # create the url for iRefIndex depending on species_id and version
  url <- sprintf("https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/%s.mitab.%s.txt.zip",
                 species_id,
                 version)
  # return the url
  return(url)
}


# MINT -----------

#' #' urlmaker_mint()
#'
#' @param species from which species does the data come from, default value = "Homo sapiens"
#' @param version version of data files in MINT, default value = "current" #CURRENTVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_mint <- urlmaker_mint(species = "Homo Sapiens")
#' url_mint

urlmaker_mint <- function (species = "Homo Sapiens", # default value = "Homo Sapiens"
                           version = "current") { # default value = current , can not fetch previous versions #CURRENTVERSION

  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1


  # create a list that contains species names and corresponding names in the url

  info_species_mint <- list( "all organisms" = "*",
                        "Homo Sapiens" = "species:human",
                        "Mus Musculus" = "species:mouse",
                        "Drosophila Melanogaster" = "species:fruit%20fly",
                        "Saccharomyces Cerevisiae" = "species:yeast")

  # check that the value for species is listed in MINT
  if (!species %in% names(info_species_mint)) {
    stop("Species not found as specified by MINT, ",
         "please check some valid entries in info_species_mint or on the website 'https://mint.bio.uniroma2.it/index.php/download/'")
    }


  # fetch the species_name from info_species_mint
  species_name <- info_species_mint[[species]]

  # create the url for MINT depending on version and species_name
  url <- sprintf( "http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/%s/search/query/%s",
                  version,
                  species_name)

  # return the url
  return(url)
}



# GeneMania -------------

#' urlmaker_genemania()
#'
#' @param species types listed in list_species_genemania depending on current version, default value = "Homo_sapiens"
#' @param version version of GeneMania , default value = "current" #UPDATEVERSION
#'
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_genemania <- urlmaker_genemania(species = "Homo_sapiens",
#'                                     version = "current")
#' url_genemania
urlmaker_genemania <- function ( species = "Homo_sapiens",
                                 version = "current") { #default value = "current"

  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1


  # create a list of all species in GeneMania

  info_species_genemania <- c("Arabidopsis_thaliana",
                              "Caenorhabditis_elegans",
                              "Danio_rerio",
                              "Drosophila_melanogaster",
                              "Escherichia_coli",
                              "Homo_sapiens",
                              "Mus_musculus",
                              "Rattus_norvegicus",
                              "Saccharomyces_cerevisiae")

  # make sure that files stored in archive have another value in url
  if (version == "current")
    archive <- ""

  else
    archive <- "archive/"

  # check that the value for species is listed in GeneMania
  if (! species %in% info_species_genemania) {
    stop("Species not found as specified by GeneMania, ",
         "please check some valid entries in info_species_genemania or on the website 'https://genemania.org/data/current/'")
  }


  # create the link depending on the version and species
  url <- sprintf("https://genemania.org/data/%s%s/%s.COMBINED/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt",
                 archive,
                 version,
                 species)

  # return the url
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
#' url_cpdb <- urlmaker_cpdb(species = "human")
#' url_cpdb
urlmaker_cpdb <- function (species = "human") { #default value = human because at the moment (05.08.2024) only human

  list_species_cpdb <- c("human", "mouse", "yeast")

  stopifnot(is.character(species)) # make sure to type in a species name as character

  if (!(species %in% list_species_cpdb)) { # if species is not in the list
    stop("Species not found as specified by CPDB,",
         "CPDB only contains data for 'human', 'mouse' and 'yeast'") # stop function and print
  }


  url <- sprintf("http://cpdb.molgen.mpg.de/download/ConsensusPathDB_%s_PPI.gz",
                 species)
  return(url)

}
# HuRi - ONLY HUMAN ----------------

#' urlmaker_huri()
#'
#' @param species default value = "human", because this database only provides human data
#' @param type different datasets , more information on "http://www.interactome-atlas.org/about/"
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#'
#' url_huri <- urlmaker_huri(species = "human",
#'                           type = "HI-union")
#' url_huri
urlmaker_huri <- function (species = "human", # default value human, because this database only provides human data
                           type = c("HI-union", "Lit-BM")) { #recommended value = "HI-union", because it contains nearly all data from HuRi



  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(type))

  if (species != "human") { # if species is not in the list
    stop("Species not found as specified by HuRi,",
         "HuRi only contains data for 'human'") # stop function and print
  }

  # the datafile "HI-union" is an aggregate of all PPIs identified in
  # HI-I-05,
  # HI-II-14
  # HuRI
  # Venkatesan-09
  # Yu-11
  # Yang-16
  # Test space screens-19

  url <- sprintf("http://www.interactome-atlas.org/data/%s.tsv",
                 type)


  return(url)

}



# TODO


# TODO
# MatrixDB - ONLY HUMAN ------------------

#' urlmaker_matrixdb()
#'
#' @param species default value = "human", because only one version and one species at MatrixDB #UPDATEVERSION
#' @param type datasets provided by MatrixDB: "CORE" = MatrixDB manually curated interaction dataset
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_matrixdb <- urlmaker_matrixdb(type = "all")
#' url_matrixdb
urlmaker_matrixdb <- function(species = "human",
                              type = c("all", "CORE")){ #UPDATEVERSION

  stopifnot(is.character(species))                  # make sure to type in a species name as character

  # check that the value for species is listed in MatrixDB
  if (species != "human") {
    stop("Species not found as specified by MatrixDB,
         MatrixDB only provide data for 'human'")
  }

  url <- sprintf("https://matrixdb.univ-lyon1.fr/downloads/matrixdb_%s.tab.zip",
                 type)

return(url)

}





# PathwayCommons (PC) - ONLY HUMAN ----------------

#' urlmaker_pathwaycommons()
#'
#' @param species default value = "Homo sapiens", because only one species
#' @param version default value = "v12" because v14 doesn't contain all datafiles #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#'
#' url_pc <- urlmaker_pathwaycommons()
#' url_pc
urlmaker_pathwaycommons<- function( species = "human", #default value = "human", because PC mostly provides data for human  #check in: c14:<unique_id> (bio processes and participants). BioPAX URIs are not to guess; instead, they should be discovered with /search or /top_pathways
                                    version = "v12") { # default value = "v12", because in "v14" not all datafiles are updated already #UPDATEVERSION

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1)  # make sure to type in a version with the length == 1

  # check that the value for species is listed in PathwayCommons
  if (species != "human") { # if species is not in the list
    stop("Species not found as specified by PathwayCommons,",
         " PathwayCommons only contains data for 'Homo sapiens'") # stop function and print
  }

  # create url depending on the version
  url <- sprintf("https://download.baderlab.org/PathwayCommons/PC2/%s/PathwayCommons12.All.hgnc.txt.gz",
                 version)

  # return the url
  return (url)

}








# Reactome - SPECIESDEFINITION ---------------
#' urlmaker_reactome()
#'
#' @param version default value = "current" #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' url_reactome <- urlmaker_reactome()
#' url_reactome

urlmaker_reactome <- function(version = "current"){ #SPECIESDEFINITION

  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  url <- "https://reactome.org/download/current/interactors/reactome.all_species.interactions.psi-mitab.txt"
  return(url)
}


# Innate DB - SPECIESDEFINITION -------------------


#' urlmaker_innatedb()
#'
#' @param url innatedb doesn't provide information about the version or the species in their url, that why the default value of url is url
#' @param version parameter set to current version 5.4 #UPDATEVERSION
#'
#' @return url
#' @export
#'
#' @examples
#' url_innatedb <- urlmaker_innatedb()
#' url_innatedb
urlmaker_innatedb <- function(url = url,
                              version = "5.4"){ # default value = "5.4" #UPDATEVERSION

  url <- "https://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz"

  return(url)
}  #SPECIESDEFINITION


# BIOGRID - SPECIESDEFINITION -------------------------------------

#' urlmaker_biogrid()
#'
#' @param version version of the data files in BioGRID
#'
#' @return url returns the corresponding url set by params #SPECIESDEFINITION
#' @export
#'
#' @examples
#' url_biogrid <- urlmaker_biogrid()
#' url_biogrid

urlmaker_biogrid <- function(version = "4.4.236") { # default value = "4.4.236" (August 2024) #SPECIESDEFINITION

  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  #create the url to download a zip file for ALL species depending on the version
  url <- sprintf("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-%s/BIOGRID-ORGANISM-%s.tab3.zip",
                 version,
                 version)

  #return the url
  return (url)

}



# CORUM - SPECIESDEFINITION ------------------------------

## #' urlmaker_corum()
## #'
## #' @param version version of the data files in CORUM, default value = "current" (August 2024 = 28.11.2022 Corum 4.1 release)
## #'
## #' @return url returns the corresponding url set by params #SPECIESDEFINITION later on
## #' @export
## #'
## #' @examples
## #' url_corum <- urlmaker_corum()
## #' url_corum
## urlmaker_corum <- function (version = "current"){ #default value set to current, but make sure to check whcih version currrent reflects (August 2024 = 28.11.2022 Corum 4.1 release) #UPDATEVERSION
## #'
## #'   stopifnot(is.character(version))                  # make sure to type in a version as character
## #'   stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1
## #'
## #'
## #'   url <- sprintf("https://mips.helmholtz-muenchen.de/corum/download/releases/%s/allComplexes.txt.zip",
## #'                  version)
## #'   return (url)
## #'
## #' }
## #'
## #' new corum :https://mips.helmholtz-muenchen.de/corum/download




# IntAct - SPECIESDEFINITION ---------------------------------

#' urlmaker_intact()
#'
#' @param version version of the data files in IntAct, default value = "current" (August 2024  = 2024-05-23 18:09	6.6G)
#'
#' @return url returns the corresponding url set by params #SPECIESDEFINITION later on
#' @export
#'
#' @examples
#' url_intact <- urlmaker_intact(version = "current")
#' url_intact

urlmaker_intact <- function(version = "current") { # default value for version, because Intact provides this file only for current (August 2024 = 2024-05-23 18:09	6.6G)

  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  # as there is only the current version of this file we can't change the link so we need
  if (!version == "current")
    stop("make sure that you use the current version of the intact.txt file, as there is only one url")

  # create url depending on version
  url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/intact/%s/psimitab/intact.txt",
                 version)

  #return the url
  return(url)

}


