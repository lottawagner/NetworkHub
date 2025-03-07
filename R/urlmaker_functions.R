

# Comment on comments
# UPDATEVERSION - for databases where version is not defined in url
# SPECIESDEFINITION - for databases where species is not defined in url
# CURRENTVERSION - for databases only providing data for current version (no archive)


#' URL maker for STRINGDB
#' 
#' Creating the URL to access the resources on [STRINGDB](https://string-db.org/) 
#' 
#' @details
#' The urlmaker function will automatically handle the conversion between human
#' readable identifiers and the identifier used internally by the selected resource
#'
#' @param type Character string, which data file do you want to download? 
#' Can be one of "PPI" or "protein_info"
#' @param species Character string, from which species does the data come from.
#' Defaults to "Homo sapiens" for human
#' @param version Character string, specifying the version of the data files in 
#' stringdb - defaults to "12.0" for the latest release (as of Feb 2025)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#'
#' @importFrom utils read.delim
#' 
#' @references https://string-db.org/
#' 
#' @family urlmakers
#'
#' @examples
#' url_stringdb <- urlmaker_stringdb(
#'   type = "PPI",
#'   species = "Homo sapiens",
#'   version = "12.0"
#' )
#'
#' url_stringdb
urlmaker_stringdb <- function(type = c("PPI", "protein_info"),
                              species = "Homo sapiens",
                              version = "12.0") {
  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  url_species_stringdb <- sprintf(
    "https://stringdb-downloads.org/download/species.v%s.txt",
    version
  )

  # read.delim the data from the species text file (columns separated using a delimiter)
  df_species <- read.delim(url(url_species_stringdb))
  info_species <- df_species

  # check that the value for species is listed in StringDB
  if (!(species %in% info_species$official_name_NCBI)) {
    stop(
      "Species not found as specified by STRINGDB, ",
      "please check some valid entries by running `info_species_stringdb()`"
    )
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
      "https://stringdb-downloads.org/download/protein.aliases.v%s/%s.protein.aliases.v%s.txt.gz",
      version,
      species_id,
      version
    )
  }
  # return value is "url"
  return(url)
}


#' URL maker for HINT
#' 
#' Creating the URL to access the resources on [HINT](https://hint.yulab.org/) 
#' (High-quality Interactomes)
#' 
#' @param type Character string, specifying which interaction types to retrieve
#' from HINT, default value = "binary"
# "binary" = binary
# "cocomp" = co-complex
# "lcb" = literature curated binary
# "lcc" = literature curated co-complex
#' @param species Character string, types listed in list_species_hint depending 
#' on current version, default value = "HomoSapiens" 
#' @param version Character string, version of HINT, always specified as 
#' `year-month` (2024-06, 2020-08, ...)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://hint.yulab.org/
#' 
#' @family urlmakers
#'
#' @examples
#' url_hint <- urlmaker_hint(
#'   type = "binary",
#'   species = "HomoSapiens",
#'   version = "2024-06"
#' )
#' 
#' url_hint
urlmaker_hint <- function(type = "binary", # default value for type = "binary"
                          species = "HomoSapiens", # default value for species #UPDATEVERSION
                          version) { # version is always written as year - month ( 2024-06, 2020-08, ...) #UPDATEVERSION

  # avoid errors
  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1


  # define the opportunities for species in HINT

  list_species_hint <- c(
    "HomoSapiens",
    "SaccharomycesCerevisiae",
    "SchizosaccharomycesPombe",
    "MusMusculus",
    "DrosophilaMelanogaster",
    "CaenorhabditisElegans",
    "ArabidopsisThaliana",
    "EscherichiaColi",
    "RattusNorvegicus",
    "OryzaSativa"
  )

  # list species is actualized for version HINT "2024-06"
  # UPDATEVERSION

  # check that the value for species is listed in HINT

  if (!(species %in% list_species_hint)) { # if species is not in the list
    stop(
      "Species not found as specified by HINT,",
      "please check some valid entries by running `list_species_hint`"
    ) # stop function and print
  }

  # define the opportunities for the four interaction types in HINT:

  # "binary" = binary
  # "cocomp" = co-complex
  # "lcb" = literature curated binary
  # "lcc" = literature curated co-complex

  type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))

  # create the url depending on the type, version and species

  url <- sprintf(
    "https://hint.yulab.org/download-raw/%s/%s_%s_hq.txt",
    version,
    species,
    type
  )

  return(url)
}


#' URL maker for FunCoup
#' 
#' Creating the URL to access the resources on [FunCoup](https://funcoup.org/)
#' 
#' @param version Character string, specifying the version of FunCoup to 
#' retrieve, default value = "6.0"
#' @param species Character string, types listed in `list_species_funcoup` 
#' depending on current version. Defaults to "H.sapiens"
#' @param type Character string, specifying the interaction types in FunCoup 
#' (for current version compact recommended, for older versions only full possible)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://funcoup.org/
#' 
#' @family urlmakers
#'
#' @examples
#' url_funcoup <- urlmaker_funcoup(
#'   version = "6.0",
#'   species = "B.taurus",
#'   type = "full"
#' )
#' 
#' url_funcoup
urlmaker_funcoup <- function(version = "6.0", # default value = "6.0"
                             species = "H.sapiens", # default value = "H.sapiens", value as type first letter (capital) of first name, full second name with small letters
                             type = c("compact", "full")) { # for current version compact recommended, for older versions only full possible

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # define the opportunities for species in FunCoup 5.0

  message("NOTE: FunCoup provides different species for the differnet version and uses different names. Check on <https://funcoup.org/archive/> for the correct definition.")

  # using sprintf to define the url of FunCoup by looking at the corresponding version and species

  # current version is stored on website in download folder #UPDATEVERSION
  if (version == "6.0") {
    list_species_funcoup_6.0 <- c(
      "A.thaliana",
      "B.subtilis",
      "B.taurus",
      "C.elegans",
      "C.familiaris",
      "C.intestinalis",
      "D.discoideum",
      "D.melanogatser",
      "D.rerio",
      "E.coli",
      "G.gallus",
      "H.sapiens",
      "M.jannaschii",
      "M.musculus",
      "M.tuberculosis",
      "O.sativa",
      "P.falciparum",
      "R.norvegicus",
      "S.cerevisiae",
      "S.pombe",
      "S.scrofa",
      "S.solfataricus",
      "SARS-CoV-2"
    )

    if (!species %in% list_species_funcoup_6.0) {
      stop(
        "Species not found as specified by Funcoup version 6.0,",
        "please check some valid entries on the webiste 'https://funcoup.org/archive/'"
      )
    }

    url <- sprintf(
      "https://funcoup.org/download/network&FC%s_%s_%s.gz",
      version,
      species,
      type
    )
  }

  if (version == "5.0") {
    list_species_funcoup_5.0 <- c(
      "A.thaliana",
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
      "S.cerevisiae",
      "S.pombe",
      "S.scrofa",
      "S.solfataricus"
    )

    if (!species %in% list_species_funcoup_5.0) {
      stop(
        "Species not found as specified by Funcoup version 5.0,",
        "please check some valid entries on the webiste 'https://funcoup.org/archive/'"
      )
    }

    url <- sprintf(
      "https://funcoup.org/download/archive&FC%s_%s_%s.gz",
      version,
      species,
      type
    )

    # archived versions are in archive folder #UPDATEVERSION
  }

  if (version == "4.1" || version == "4.0") {
    list_species_funcoup_4._ <- c(
      "A.thaliana",
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
      "S.cerevisiae",
      "S.pombe"
    )
    if (!species %in% list_species_funcoup_4._) {
      stop(
        "Species not found as specified by Funcoup version 4.0 or 4.1,",
        "please check some valid entries on the webiste 'https://funcoup.org/archive/'"
      )
    }

    url <- sprintf(
      "https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_full.gz",
      version,
      version,
      species
    )

    if (type == "compact") {
      message("NOTE: For version 4.0 & 4.1: Only the full datafile is provided by funcoup.")
    }
  }

  if (version == "3.0") {
    list_species_funcoup_3.0 <- c(
      "A.thaliana",
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
      "S.cerevisiae"
    )

    if (!species %in% list_species_funcoup_3.0) {
      stop(
        "Species not found as specified by Funcoup version 3.0,",
        "please check some valid entries on the webiste 'https://funcoup.org/archive/'"
      )
    }

    url <- sprintf(
      "https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=FC%s_%s_%s.gz",
      version,
      version,
      species,
      type
    )
  }

  if (version == "2.0") {
    list_species_funcoup_2.0 <- c(
      "athaliana",
      "celegans",
      "cfamiliaris",
      "cintestinalis",
      "dmelanogaster",
      "drerio",
      "ggallus",
      "hsapiens",
      "mmusculus",
      "rnorvegicus",
      "scerevisiae"
    )

    if (!species %in% list_species_funcoup_2.0) {
      stop(
        "Species not found as specified by Funcoup version 2.0,",
        "please check some valid entries on the webiste 'https://funcoup.org/archive/'"
      )
    }

    url <- sprintf(
      "https://funcoup.org/archive/download.action?type=archive&instanceID=24480085&version=FunCoup-%s&fileName=%s.%s.pfc01.tsv.gz",
      version,
      species,
      type
    )
  }

  if (version == "1.0") {
    stop("No url provided by NetworkHub for version 1.0")
  }

  return(url)
}




#' URL maker for IID
#' 
#' Creating the URL to access the resources on [IID](https://iid.ophid.utoronto.ca/)
#' (Integrated Interactions Database)
#' 
#' @param species Character string, with species listed in `list_species_iid` 
#' depending on current version, default value = "human"
#' @param version Character string, version of IID , default value = "2021-05"
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://iid.ophid.utoronto.ca/
#'
#' @family urlmakers
#' 
#' @examples
#' url_iid <- urlmaker_iid(species = "mouse")
#' 
#' url_iid
urlmaker_iid <- function(species = "human", #
                         version = "2021-05") { # version of IID not updated since

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1


  # define the opportunities for species in IID

  list_species_iid <- c(
    "alpaca",
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
    "yeast"
  )
  # list species is actualized for version IID version 2021-05
  # UPDATEVERSION

  # check that the value for species is listed in IID
  if (!species %in% list_species_iid) {
    stop(
      "Species not found as specified by IID,",
      "please check some valid entries of `list_species_iid` and on the webiste 'https://iid.ophid.utoronto.ca/search_by_proteins/'"
    )
  }

  # create the url depending on the species
  url <- sprintf(
    "https://iid.ophid.utoronto.ca/static/download/%s_annotated_PPIs.txt.gz",
    species
  )

  return(url)
}


#' URL maker for iRefIndex
#' 
#' Creating the URL to access the resources on [iRefIndex](https://irefindex.vib.be/)
#' (Integrated Interactions Database)
#' 
#' @param species Character string, from which species does the data come from, 
#' default value = "Homo sapiens"
#' @param version Character string, version of data files in iRefIndex, 
#' default value = "08-28-2023" #UPDATEVERSION
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://irefindex.vib.be/
#'
#' @family urlmakers
#' 
#' @examples
#' url_irefindex <- urlmaker_irefindex(
#'   species = "Homo sapiens",
#'   version = "08-28-2023"
#' )
#' url_irefindex
urlmaker_irefindex <- function(species,
                               version = "08-28-2023") {

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # create a list that contains species names and corresponding IDs

  list_species_irefindex <- c(
    "Homo sapiens",
    "Mus musculus",
    "Saccharomyces cerevisiae S288C",
    "Escherichia",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Drosophila melanogaster",
    "Caenorhabditis elegans"
  )

  species_id_irefindex <- c(
    "9606",
    "10090",
    "559292",
    "562",
    "10116",
    "4932",
    "7227",
    "6239"
  )


  irefindex_db_annotations <- data.frame(
    species_irefindex = list_species_irefindex,
    species_id = species_id_irefindex,
    row.names = list_species_irefindex
  )

  # check that the value for species is listed in iRefIndex
  if (!species %in% list_species_irefindex) {
    stop(
      "Species not found as specified by iRefIndex, ",
      "please check some valid entries in info_species_irefindex_id or on the website 'https://irefindex.vib.be/wiki/index.php/README_MITAB2.6_for_iRefIndex_20.0#Column_number:_29_.28Host_organism_taxid.29`"
    )
  }

  # fetch the species_id from info_species_irefindex
  species_row <- irefindex_db_annotations[irefindex_db_annotations$species_irefindex == species, ]
  species_id <- species_row$species_id

  # create the url for iRefIndex depending on species_id and version
  url <- sprintf(
    "https://storage.googleapis.com/irefindex-data/archive/release_20.0/psi_mitab/MITAB2.6/%s.mitab.%s.txt.zip",
    species_id,
    version
  )

  return(url)
}


#' URL maker for MINT
#' 
#' Creating the URL to access the resources on [MINT](https://mint.bio.uniroma2.it/)
#' (The Molecular INTeraction Database)
#'
#' @param species Character string, from which species does the data come from, 
#' default value = "Homo Sapiens"
#' @param version Character string, version of data files in MINT, defaults to 
#' "current" (as it is specified by the MINT curators)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://mint.bio.uniroma2.it/
#'
#' @family urlmakers
#'
#' @examples
#' url_mint <- urlmaker_mint(species = "Homo Sapiens")
#' url_mint
urlmaker_mint <- function(species = "Homo Sapiens", # default value = "Homo Sapiens"
                          version = "current") { # default value = current , can not fetch previous versions #CURRENTVERSION

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1


  # create a list that contains species names and corresponding names in the url

  info_species_mint <- list(
    "all organisms" = "*",
    "Homo Sapiens" = "species:human",
    "Mus Musculus" = "species:mouse",
    "Drosophila Melanogaster" = "species:fruit%20fly",
    "Saccharomyces Cerevisiae" = "species:yeast"
  )

  # check that the value for species is listed in MINT
  if (!species %in% names(info_species_mint)) {
    stop(
      "Species not found as specified by MINT, ",
      "please check some valid entries in info_species_mint or on the website 'https://mint.bio.uniroma2.it/index.php/download/'"
    )
  }

  # fetch the species_name from info_species_mint
  species_name <- info_species_mint[[species]]

  # create the url for MINT depending on version and species_name
  url <- sprintf(
    "http://www.ebi.ac.uk/Tools/webservices/psicquic/mint/webservices/%s/search/query/%s",
    version,
    species_name
  )

  return(url)
}


#' URL maker for GeneMania
#' 
#' Creating the URL to access the resources on [GeneMania](https://genemania.org/)
#'
#' @param species Character string, types listed in `list_species_genemania` 
#' depending on current version, default value = "Homo_sapiens"
#' @param version Character string, version of GeneMania , defaults to "current"
#'
#'
#' @return A character string with the URL as requested
#' 
#' @references https://genemania.org/
#' 
#' @export
#'
#' @family urlmakers
#' 
#' @examples
#' url_genemania <- urlmaker_genemania(
#'   species = "Homo_sapiens",
#'   version = "current"
#' )
#' 
#' url_genemania
urlmaker_genemania <- function(species = "Homo_sapiens",
                               version = "current") { # default value = "current"

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # create a list of all species in GeneMania
  info_species_genemania <- c(
    "Arabidopsis_thaliana",
    "Caenorhabditis_elegans",
    "Danio_rerio",
    "Drosophila_melanogaster",
    "Escherichia_coli",
    "Homo_sapiens",
    "Mus_musculus",
    "Rattus_norvegicus",
    "Saccharomyces_cerevisiae"
  )

  # make sure that files stored in archive have another value in url
  if (version == "current") {
    archive <- ""
  } else {
    archive <- "archive/"
  }

  # check that the value for species is listed in GeneMania
  if (!species %in% info_species_genemania) {
    stop(
      "Species not found as specified by GeneMania, ",
      "please check some valid entries in info_species_genemania or on the website 'https://genemania.org/data/current/'"
    )
  }


  # create the link depending on the version and species
  url <- sprintf(
    "https://genemania.org/data/%s%s/%s.COMBINED/COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt",
    archive,
    version,
    species
  )

  return(url)
}



#' URL maker for HuRI
#' 
#' Creating the URL to access the resources on 
#' [HuRI](http://www.interactome-atlas.org/) (The Human Reference Protein 
#' Interactome)
#'
#' @param species Character string, default value = "human", because this 
#' database only provides human data
#' @param type Character string, different datasets , more information on 
#' "http://www.interactome-atlas.org/about/" - defaults to "HI-union"
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references http://www.interactome-atlas.org/
#' 
#' @family urlmakers
#'
#' @examples
#'
#' url_huri <- urlmaker_huri(
#'   species = "human",
#'   type = "HI-union"
#' )
#' url_huri
urlmaker_huri <- function(species = "human", 
                          type = c("HI-union", "Lit-BM") # recommended value = "HI-union", because it contains nearly all data from HuRI
                          ) { 

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(type))

  if (species != "human") { # if species is not in the list
    stop(
      "Species not found as specified by HuRi,",
      "HuRi only contains data for 'human'"
    ) # stop function and print
  }

  # the datafile "HI-union" is an aggregate of all PPIs identified in
  # HI-I-05,
  # HI-II-14
  # HuRI
  # Venkatesan-09
  # Yu-11
  # Yang-16
  # Test space screens-19

  url <- sprintf(
    "http://www.interactome-atlas.org/data/%s.tsv",
    type
  )

  return(url)
}



#' URL maker for MatrixDB
#' 
#' Creating the URL to access the resources on 
#' [MatrixDB](https://matrixdb.univ-lyon1.fr/)
#'
#' @param species Character string, default value = "human", because only one 
#' version and one species at MatrixDB #UPDATEVERSION
#' @param type Character string, datasets provided by MatrixDB: "CORE" 
#' = MatrixDB manually curated interaction dataset. Defaults to "all"
#' @param version Character string, specifying the version number. Defaults to 
#' the recently updated "4.0" (as this was recently added to the URL)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://matrixdb.univ-lyon1.fr/
#' 
#' @family urlmakers
#'
#' @examples
#' url_matrixdb <- urlmaker_matrixdb(type = "CORE")
#'
#' url_matrixdb
urlmaker_matrixdb <- function(species = "human",
                              type = c("all", "CORE"),
                              version = "4_0") { # UPDATEVERSION

  stopifnot(is.character(species)) # make sure to type in a species name as character

  # check that the value for species is listed in MatrixDB
  if (species != "human") {
    stop("Species not found as specified by MatrixDB,
         MatrixDB only provide data for 'human'")
  }

  url <- sprintf(
    "https://matrixdb.univ-lyon1.fr/downloads/matrixdb_%s_%s.tab.zip",
    type,
    version
  )

  return(url)
}


#' URL maker for PathwayCommons
#' 
#' Creating the URL to access the resources on 
#' [PathwayCommons](https://www.pathwaycommons.org/)
#'
#' @param species Character string, default value = "Homo sapiens", because 
#' only one species is available
#' @param version Character string, default value = "v12" (v14 doesn't contain 
#' all datafiles as of Feb 2025)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://www.pathwaycommons.org/
#' 
#' @family urlmakers
#'
#' @examples
#'
#' url_pc <- urlmaker_pathwaycommons()
#' 
#' url_pc
urlmaker_pathwaycommons <- function(species = "human", # default value = "human", because PC mostly provides data for human  #check in: c14:<unique_id> (bio processes and participants). BioPAX URIs are not to guess; instead, they should be discovered with /search or /top_pathways
                                    version = "v12") { # default value = "v12", because in "v14" not all datafiles are updated already #UPDATEVERSION

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # check that the value for species is listed in PathwayCommons
  if (species != "human") { # if species is not in the list
    stop(
      "Species not found as specified by PathwayCommons,",
      " PathwayCommons only contains data for 'Homo sapiens'"
    ) # stop function and print
  }

  # create url depending on the version
  url <- sprintf(
    "https://download.baderlab.org/PathwayCommons/PC2/%s/PathwayCommons12.All.hgnc.txt.gz",
    version
  )

  return(url)
}


#' URL maker for HIPPIE
#' 
#' Creating the URL to access the resources on 
#' [HIPPIE](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/)
#'
#' @param species Character string, default value = "Homo_sapiens", because 
#' this database only provides human data
#' @param version Character string, default value = "current", 
#' version of the database ... #UPDATEVERSION
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/
#' 
#' @family urlmakers
#'
#' @examples
#'
#' url_hippie <- urlmaker_hippie(
#'   species = "Homo_sapiens",
#'   version = "current"
#' )
#' url_hippie
#'
urlmaker_hippie <- function(species = "Homo_sapiens", # default value human, because this database only provides human data
                            version = "current") { # default value current

  stopifnot(is.character(species)) # make sure to type in a species name as character
  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  if (species != "Homo_sapiens") { # if species is not in the list
    stop(
      "Species not found as specified by HIPPIE,",
      "HIPPIE only contains data for 'Homo_sapiens'"
    ) # stop function and print
  }

  url <- sprintf(
    "https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/HIPPIE-%s.mitab.txt",
    version
  )

  return(url)
}


#' URL maker for Reactome
#' 
#' Creating the URL to access the resources on [Reactome](https://reactome.org/)
#'
#' @param version Character string, default value = "current" #UPDATEVERSION
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://reactome.org/
#' 
#' @family urlmakers
#'
#' @examples
#' url_reactome <- urlmaker_reactome()
#' 
#' url_reactome
urlmaker_reactome <- function(version = "current") { # SPECIESDEFINITION

  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  url <- "https://reactome.org/download/current/interactors/reactome.all_species.interactions.psi-mitab.txt"
  
  return(url)
}


#' URL maker for InnateDB
#' 
#' Creating the URL to access the resources on 
#' [InnateDB](https://www.innatedb.com/)
#'
#' @param url Character string, innatedb doesn't provide information about the version or the species in their url, that why the default value of url is url
#' @param version Character string, parameter set to current version 5.4 #UPDATEVERSION
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://www.innatedb.com/
#'
#' @family urlmakers
#' 
#' @examples
#' url_innatedb <- urlmaker_innatedb()
#' 
#' url_innatedb
urlmaker_innatedb <- function(url = url,
                              version = "5.4") { # default value = "5.4" #UPDATEVERSION

  url <- "https://www.innatedb.com/download/interactions/innatedb_ppi.mitab.gz"

  return(url)
}


#' URL maker for BioGRID
#' 
#' Creating the URL to access the resources on [BioGRID](https://thebiogrid.org/)
#'
#' @param version version of the data files in BioGRID
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://thebiogrid.org/
#' 
#' @family urlmakers
#'
#' @examples
#' url_biogrid <- urlmaker_biogrid()
#' 
#' url_biogrid
urlmaker_biogrid <- function(version = "4.4.238") { # default value = "4.4.238" (October 2024) #SPECIESDEFINITION

  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # create the url to download a zip file for ALL species depending on the version
  url <- sprintf(
    "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-%s/BIOGRID-ALL-%s.tab3.zip",
    version,
    version
  )

  return(url)
}


#' URL maker for IntAct
#' 
#' Creating the URL to access the resources on 
#' [IntAct](https://www.ebi.ac.uk/intact/)
#'
#' @param version Character string, version of the data files in IntAct, 
#' default value = "current" (August 2024  = 2024-05-23 18:09	6.6G)
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references https://www.ebi.ac.uk/intact/
#' 
#' @family urlmakers
#'
#' @examples
#' url_intact <- urlmaker_intact(version = "current")
#' 
#' url_intact
urlmaker_intact <- function(version = "current") { # default value for version, because Intact provides this file only for current (August 2024 = 2024-05-23 18:09	6.6G)

  stopifnot(is.character(version)) # make sure to type in a version as character
  stopifnot(length(version) == 1) # make sure to type in a version with the length == 1

  # as there is only the current version of this file we can't change the link so we need
  if (!version == "current") {
    stop("make sure that you use the current version of the intact.txt file, as there is only one url")
  }

  # create url depending on version
  url <- sprintf(
    "https://ftp.ebi.ac.uk/pub/databases/intact/%s/psimitab/intact.zip",
    version
  )

  return(url)
}


#' URL maker for ConsensusPathDB
#' 
#' Creating the URL to access the resources on 
#' [ConsensusPathDB](http://cpdb.molgen.mpg.de/)
#'
#' @param species Character string, specifying the species of interesting,
#' defaulting to  "human". 
#'
#' @return A character string with the URL as requested
#' 
#' @export
#' 
#' @references http://cpdb.molgen.mpg.de/
#' 
#' @family urlmakers
#'
#' @examples
#' url_cpdb <- urlmaker_cpdb(species = "human")
#' 
#' url_cpdb
urlmaker_cpdb <- function(species = "human") { # default value = human because at the moment (05.08.2024) only human
  
  list_species_cpdb <- c("human", "mouse", "yeast")
  
  stopifnot(is.character(species)) # make sure to type in a species name as character
  
  if (!(species %in% list_species_cpdb)) { # if species is not in the list
    stop(
      "Species not found as specified by CPDB,",
      "CPDB only contains data for 'human', 'mouse' and 'yeast'"
    ) # stop function and print
  }
  
  url <- sprintf(
    "http://cpdb.molgen.mpg.de/download/ConsensusPathDB_%s_PPI.gz",
    species
  )
  
  return(url)
}
