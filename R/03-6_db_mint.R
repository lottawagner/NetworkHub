# get_networkdata_mint() -----------

#' get_networkdata_mint()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in mint
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation default value set to TRUE (MINT annotation dataframe contains Uniprot, GeneSymbol and Ensembl)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_mint
#'
#' @importFrom vroom vroom
#' @importFrom stringr str_extract
#'
#' @export
#'
#' @examples
#' \donttest{
#' db_mint_df <- get_networkdata_mint(
#'   species = "Homo Sapiens",
#'   version = "current"
#' )
#'
#' db_mint_df
#' }

get_networkdata_mint <- function(species,
                                 version = "current",
                                 cache = TRUE,
                                 add_annotation = TRUE,
                                 ...) {



  # list species is actualized for version mint "current"
  # UPDATEVERSION

  # check that the value for species is listed in mint

  if (!(species %in% list_species_mint)) { # if species is not in the list
    stop("Species not found as specified by mint,",
         "please check some valid entries by running `list_species_mint`") # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "mint_",
    species,
    "_v",
    version
  ) # definition of the resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" mint url
    mint_url <-
      urlmaker_mint(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_mint

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = mint_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_mint <- vroom::vroom(network_file, col_names = FALSE)
  #ppis_mint <- head(read.delim(network_file, sep = " "))
  #
  message(dim(ppis_mint))

  colnames(ppis_mint) <- c("Identifier_A", "Identifier_B", "Alternative identifier for interactor A", "Alternative identifier for interactor B", "Aliases for A", "Aliases for B", "method", "First author", "pubmed", "NCBI Taxonomy identifier for interactor A","NCBI Taxonomy identifier for interactor B", "interaction_type", "Source databases and identifiers", "Interaction identifier(s) in the corresponding source database", "Confidence score")

  if (add_annotation) {

    ppi_mint_df_annotated <- ppis_mint

    #UniProt
    ppi_mint_df_annotated$Uniprot_A <- str_extract(ppi_mint_df_annotated$Identifier_A, "uniprotkb:([A-Z0-9]+)")
    ppi_mint_df_annotated$Uniprot_A <- gsub("uniprotkb:", "", ppi_mint_df_annotated$Identifier_A)

    ppi_mint_df_annotated$Uniprot_B <- str_extract(ppi_mint_df_annotated$Identifier_B, "uniprotkb:([A-Z0-9]+)")
    ppi_mint_df_annotated$Uniprot_B <- gsub("uniprotkb:", "", ppi_mint_df_annotated$Identifier_B)


    #GeneSymbol
    ppi_mint_df_annotated$GeneSymbol_A <- str_extract(ppi_mint_df_annotated$`Aliases for A`, "uniprotkb:([^\\(]+)\\(gene name\\)")
    ppi_mint_df_annotated$GeneSymbol_A <- gsub("uniprotkb:", "", ppi_mint_df_annotated$GeneSymbol_A)
    ppi_mint_df_annotated$GeneSymbol_A <- gsub("\\(gene name\\)", "", ppi_mint_df_annotated$GeneSymbol_A)

    ppi_mint_df_annotated$GeneSymbol_B <- str_extract(ppi_mint_df_annotated$`Aliases for B`, "uniprotkb:([^\\(]+)\\(gene name\\)")
    ppi_mint_df_annotated$GeneSymbol_B <- gsub("uniprotkb:", "", ppi_mint_df_annotated$GeneSymbol_B)
    ppi_mint_df_annotated$GeneSymbol_B <- gsub("\\(gene name\\)", "", ppi_mint_df_annotated$GeneSymbol_B)

    #Ensembl
    ppi_mint_df_annotated$Ensembl_A <- strsplit(ppi_mint_df_annotated$`Alternative identifier for interactor A`, "\\|") # split entry by "|" -> you get a list for each row
    ppi_mint_df_annotated$Ensembl_A <- lapply(ppi_mint_df_annotated$Ensembl_A, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
      ensembl_entry <- grep("^ensembl:", x, value = TRUE) #search for the entry (x) in every list that starts (^) with "ensembl"
      ensembl_entry <- gsub("ensembl:", "", ensembl_entry) #remove prefix "ensembl" by using gsub()
      return(ensembl_entry)
    })

    ppi_mint_df_annotated$Ensembl_B <- strsplit(ppi_mint_df_annotated$`Alternative identifier for interactor B`, "\\|") # split entry by "|" -> you get a list for each row
    ppi_mint_df_annotated$Ensembl_B <- lapply(ppi_mint_df_annotated$Ensembl_B, function(x) { # with lapply you can iterate through every entry in a list by using an anonymous function(x), x is each element of a list
      ensembl_entry <- grep("^ensembl:", x, value = TRUE) #search for the entry (x) in every list that starts (^) with "ensembl"
      ensembl_entry <- gsub("ensembl:", "", ensembl_entry) #remove prefix "ensembl" by using gsub()
      return(ensembl_entry)
    })

    #Entrez - MINT doesn't provide info about entrez number

    return(ppi_mint_df_annotated)
  }

  if (!add_annotation) {
    return(ppis_mint)
  }
}


# outside of function ----------

list_species_mint <- c ("Homo Sapiens",
                        "Mus Musculus",
                        "Drosophila Melanogaster",
                        "Saccharomyces Cerevisiae")










