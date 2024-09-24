
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
    url <- sprintf("https://funcoup.org/downloads/download.action?type=network&instanceID=24480085&fileName=FC%s_%s_compact.gz",
                 version,
                 species)

  # archived versions are in archive folder #UPDATEVERSION
  } else if (version != "5.0") {
    url <- sprintf("https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_compact.gz",
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
#' urlmaker_irefindex <- function(species = "Homo sapiens",
#'                                 version = "2023-08-29")
#'
urlmaker_irefindex <- function(species,
                               version = "2023-08-29"){ #default value for version = "2023-08-29" #UPDATEVERSION

  stopifnot(is.character(species))                  # make sure to type in a species name as character
  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  # create a list that contains species names and corresponding IDs
  info_species_irefindex <- list(
    "Homo sapiens" = "9606",
    "Mus musculus" = "10090",
    "Saccharomyces cerevisiae" = "559292",
    "Escherichia" = "562",
    "Rattus norvegicus" = "10116",
    "Saccharomyces cerevisiae" = "4932",
    "Drosophila melanogaster" = "7227",
    "Caenorhabditis elegans" = "6239"
  )

  # check that the value for species is listed in iRefIndex
  if (!species %in% names(info_species_irefindex)) {
    stop("Species not found as specified by iRefIndex, ",
         "please check some valid entries in info_species_irefindex or on the website 'https://irefindex.vib.be/wiki/index.php/README_MITAB2.6_for_iRefIndex_20.0#Column_number:_29_.28Host_organism_taxid.29`")
  }

  # fetch the species_id from info_species_irefindex
  species_id <- info_species_irefindex[[species]]

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
#' urlmaker_mint(species = "Homo Sapiens")
#'

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
#' urlmaker_genemania( species = "Homo_sapiens",
#'                      version = "current" )
#'
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




# Reactome ---------------
#' urlmaker_reactome()
#'
#' @param version default value = "current" #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_reactome()

urlmaker_reactome <- function(version = "current"){ #SPECIESDEFINITION

  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1

  url <- "https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt"
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

  if (!(species %in% list_species_cpdb)) { # if species is not in the list
    stop("Species not found as specified by CPDB,",
         "CPDB only contains data for 'human', 'mouse' and 'yeast'") # stop function and print
  }


  url <- sprintf("http://cpdb.molgen.mpg.de/download/ConsensusPathDB_%s_PPI.gz",
                 species)
  return(url)

}
# HuRi - ONLY HUMAN ----------------

#' urlmaker_huri
#'
#' @param species default value = "human", because this database only provides human data
#' @param type different datasets , more information on "http://www.interactome-atlas.org/about/"
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#'
#' urlmaker_huri(species = "human",
#'              type = "HI-union")
urlmaker_huri <- function (species = "human", # default value human, because this database only provides human data
                           type = "HI-union") { #default value = "HI-union", because it contains nearly all data from HuRi



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

  type <- c("HI-union", "Lit-BM")

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
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#' urlmaker_matrixdb()
#'
urlmaker_matrixdb <- function(species = "human"){ #UPDATEVERSION

  stopifnot(is.character(species))                  # make sure to type in a species name as character

  # check that the value for species is listed in MatrixDB
  if (species != "human") {
    stop("Species not found as specified by MatrixDB,
         MatrixDB only provide data for 'human'")
  }

  url <- "http://matrixdb.univ-lyon1.fr/download/matrixdb_FULL.tab.gz"

return(url)

}





# PathwayCommons (PC) - ONLY HUMAN ----------------

#' urlmaker_pc()
#'
#' @param species default value = "Homo sapiens", because only one species
#' @param version default value = "v12" because v14 doesn't contain all datafiles #UPDATEVERSION
#'
#' @return url returns the corresponding url set by params
#' @export
#'
#' @examples
#'
#' urlmaker_pc()
#'
urlmaker_pc <- function( species = "Homo sapiens", #default value = "Homo sapiens", because PC mostly provides data for human  #check in: c14:<unique_id> (bio processes and participants). BioPAX URIs are not to guess; instead, they should be discovered with /search or /top_pathways
                         version = "v12") { # default value = "v12", because in "v14" not all datafiles are updated already #UPDATEVERSION

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1)  # make sure to type in a version with the length == 1

  # check that the value for species is listed in PathwayCommons
  if (species != "Homo sapiens") { # if species is not in the list
    stop("Species not found as specified by PathwayCommons,",
         " PathwayCommons only contains data for 'Homo sapiens'") # stop function and print
  }

  # create url depending on the version
  url <- sprintf("https://download.baderlab.org/PathwayCommons/PC2/%s/PathwayCommons12.All.hgnc.txt.gz",
                 version)

  # return the url
  return (url)

}








# Innate DB - SPECIESDEFINITION -------------------


urlmaker_innatedb <- function(url = url,
                              version = "5.4"){ # default value = "5.4" #UPDATEVERSION

  url <- "https://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz"
}  #SPECIESDEFINITION


# species NCBI listed in ncbi_taxid_host_organism #TODO

info_species_innatedb <- list("Mus musculus" = "10090",
                              "Homo sapiens" = "9606",
                              "Saccharomyces cerevisiae" ="4932",
                              "taxid:0",
                              "Canis lupus familiaris" = "9615",
                              "Spodoptera frugiperda" = "7108",
                              "Rattus norvegicus" = "10116",
                              "Chlorocebus aethiops" = "9534",
                              "Cricetulus griseus" = "10029",
                              "Gallus gallus" = "9031",
                              "Platyrrhini" = "9479",
                              "Oryctolagus cuniculus" = "9986",
                              "Cricetinae" = "10026",
                              "Coturnix japonica" = "93934",
                              "Bos taurus" = "9913",
                              "Escherichia phage T7" = "10760",
                              "Yeast two-hybrid vector pC-ACT.2" = "111296",
                              "unidentified baculovirus" = "10469",
                              "Chlorocebus aethiops aethiops" = "101841",
                              "Escherichia coli"  ="562",
                              "Neogale vison" = "452646",
                              "Drosophila melanogaster" = "7227",
                              "Hylobates lar" = "9580",
                              "Spodoptera frugiperda multiple nucleopolyhedrovirus" = "10455",
                              "Mesocricetus auratus" = "10036",
                              "taxid:-1",
                              "Cercopithecidae" = "9527",
                              "Mustela lutreola" =  "9666",
                              "Lepidoptera" = "7088")

# BIOGRID - SPECIESDEFINITION -------------------------------------

#' urlmaker_biogrid()
#'
#' @param version version of the data files in BioGRID
#'
#' @return url returns the corresponding url set by params #SPECIESDEFINITION
#' @export
#'
#' @examples
#' urlmaker_biogrid()

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

#' urlmaker_corum()
#'
#' @param version version of the data files in CORUM, default value = "current" (August 2024 = 28.11.2022 Corum 4.1 release)
#'
#' @return url returns the corresponding url set by params #SPECIESDEFINITION later on
#' @export
#'
#' @examples
#' urlmaker_corum()
urlmaker_corum <- function (version = "current"){ #default value set to current, but make sure to check whcih version currrent reflects (August 2024 = 28.11.2022 Corum 4.1 release) #UPDATEVERSION

  stopifnot(is.character(version))                  # make sure to type in a version as character
  stopifnot(length(version) == 1)                   # make sure to type in a version with the length == 1


  url <- sprintf("https://mips.helmholtz-muenchen.de/corum/download/releases/%s/allComplexes.txt.zip",
                 version)
  return (url)

}




# IntAct - SPECIESDEFINITION ---------------------------------

#' urlmaker_intact()
#'
#' @param version version of the data files in IntAct, default value = "current" (August 2024  = 2024-05-23 18:09	6.6G)
#'
#' @return url returns the corresponding url set by params #SPECIESDEFINITION later on
#' @export
#'
#' @examples
#' urlmaker_intact(version = "current")

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


