

# get_networkdata_reactome() -----------

#' get_networkdata_reactome()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in reactome
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with four columns (Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_reactome
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#'
#' db_reactome_df <- get_networkdata_reactome(
#'   species = "human",
#'   version = "2021-05"
#' )
#'
#' db_reactome_df
#'

get_networkdata_reactome <- function(species,
                                     version = "current",
                                     cache = TRUE,
                                     add_annotation = TRUE,
                                     ...) {

  # list species is actualized for version reactome "current" (sept 2024)
  # UPDATEVERSION

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "reactome_v",
    version,
    "_",
    species
    )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" reactome url
    reactome_url <-
      urlmaker_reactome(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_reactome

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = reactome_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_reactome <- vroom::vroom(network_file)
  #ppis_reactome <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_reactome))

  # rename columns

  #Uniprot
  colnames(ppis_reactome)[colnames(ppis_reactome) == "#ID(s) interactor A"] <- "Uniprot_A"
  colnames(ppis_reactome)[colnames(ppis_reactome) == "ID(s) interactor B"] <- "Uniprot_B"

  ppis_reactome_filtered <- ppis_reactome[(ppis_reactome$`Taxid interactor A` == species | ppis_reactome$`Taxid interactor A` == "-") &
                                          (ppis_reactome$`Taxid interactor B` == species | ppis_reactome$`Taxid interactor B` == "-"),
                                          ]

  if (add_annotation) {


    url_reactome_annotation_info <- "https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt"

    rname <- paste0(
      "reactome_",
      version,
      "_",
      species,
      "_annotation_info"
      )

    if (cache) {
      # tries to fetch from the cache
      message("Trying to fetch from cache...")
      network_file_annotation <- fetch_NetworkHub(rname)
    }

    if (!cache | is.null(network_file_annotation)) {
      # retrieves the file for the first time
      message("Downloading to cache...")
      network_file_annotation <- cache_NetworkHub(
        rname = rname,
        fpath = url_reactome_annotation_info
      )
    }

    reactome_annotation_info <- vroom::vroom(network_file_annotation)

    #rename columns Uniprot
    colnames(reactome_annotation_info)[colnames(reactome_annotation_info) == "# Interactor 1 uniprot id"] <- "Uniprot_A"
    colnames(reactome_annotation_info)[colnames(reactome_annotation_info) == "Interactor 2 uniprot id"] <- "Uniprot_B"

    #merge dataframes into new dataframe
    ppis_reactome_filtered_annotated <- merge(ppis_reactome_filtered, reactome_annotation_info, by = c("Uniprot_A", "Uniprot_B"), all.x = TRUE)
    ppis_reactome_filtered_annotated <- ppis_reactome_filtered_annotated[
      !(ppis_reactome_filtered_annotated$`Taxid interactor A` == "-" &
          ppis_reactome_filtered_annotated$`Taxid interactor B` == "-"),
    ]

    #rename columns Entrez
    colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 1 Entrez Gene id"] <- "Entrez_A"
    colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 2 Entrez Gene id"] <- "Entrez_B"

    #rename columns Ensembl
    colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 1 Ensembl gene id"] <- "Ensembl_A"
    colnames(ppis_reactome_filtered_annotated)[colnames(ppis_reactome_filtered_annotated) == "Interactor 2 Ensembl gene id"] <- "Ensembl_B"

    #extract UniProt_id
    ppis_reactome_filtered_annotated$Uniprot_A <- str_extract(ppis_reactome_filtered_annotated$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
    ppis_reactome_filtered_annotated$Uniprot_A <- gsub("uniprotkb:", "", ppis_reactome_filtered_annotated$Uniprot_A)

    ppis_reactome_filtered_annotated$Uniprot_B <- str_extract(ppis_reactome_filtered_annotated$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
    ppis_reactome_filtered_annotated$Uniprot_B <- gsub("uniprotkb:", "", ppis_reactome_filtered_annotated$Uniprot_B)


    return(ppis_reactome_filtered_annotated)
  }

  if (!add_annotation) {
    return(ppis_reactome_filtered)
  }
}







