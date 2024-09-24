

#' # get_networkdata_iid() -----------
#'
#' #' get_networkdata_iid()
#' #'
#' #' @param species  from which species does the data come from
#' #' @param version version of the data files in iid
#' #' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' #' @param ... 	further arguments passed to or from other methods
#' #'
#' #' @return ppis_iid
#' #'
#' #' @importFrom vroom vroom
#' #' @export
#' #'
#' #' @examples
#' #'
#' #' db_iid_df <- get_networkdata_iid(
#' #'   species = "human",
#' #'   version = "2021-05"
#' #' )
#' #'
#' #' db_iid_df
#' #'
#'
#' get_networkdata_iid <- function(species,
#'                                  version,
#'                                  cache = TRUE
#'                                  ...) {
#'
#'
#'
#'   # list species is actualized for version iid "2024-06"
#'   # UPDATEVERSION
#'
#'   # check that the value for species is listed in iid
#'
#'   if (!(species %in% list_species_iid)) { # if species is not in the list
#'     stop("Species not found as specified by iid,",
#'          "please check some valid entries by running `list_species_iid`") # stop function and print
#'   }
#'
#'
#'   # type <- match.arg(type, c("binary", "cocomp", "lcb", "lcc"))
#'
#'   # buildup of the resource location for the version and all
#'   ## elegantly done in another smaller utility function
#'
#'   rname <- paste0(
#'     "iid_",
#'     species,
#'     "_v",
#'     version,
#'     "_type",
#'     type
#'   ) # definition of the resource name
#'
#'   if (cache) {
#'     # tries to fetch from the cache
#'     message("Trying to fetch from cache...")
#'     network_file <- fetch_NetworkHub(rname)
#'   } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file
#'
#'   if (!cache | is.null(network_file)) {
#'     # retrieves the file for the first time
#'     message("Downloading to cache...")
#'     # buildup from "base" iid url
#'     iid_url <-
#'       urlmaker_iid(
#'         type = type,
#'         species = species,
#'         version = version
#'       ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_iid
#'
#'     # and cache_NetworkHub to cache the file from the url source
#'     network_file <- cache_NetworkHub(
#'       rname = rname,
#'       fpath = iid_url
#'     )
#'   }
#'
#'   # read in the resource, whether cached or freshly downloaded
#'   ppis_iid <- vroom::vroom(network_file)
#'   #ppis_iid <- head(read.delim(network_file, sep = " "))
#'   #
#'   message(dim(ppis_iid))
#'
#'   colnames(ppis_iid)[colnames(ppis_iid) == "Gene_A"] <- "GeneSymbol_A"
#'   colnames(ppis_iid)[colnames(ppis_iid) == "Gene_B"] <- "GeneSymbol_B"
#'
#'   if (add_annotation) {
#'     ppi_iid_df_annotated <- annotation_iid(ppi_iid = ppis_iid,
#'                                              species = species,
#'                                              version = version,
#'                                              type = type)
#'     return(ppi_iid_df_annotated)
#'   }
#'
#'
#'
#'   if (!add_annotation) {
#'     return(ppis_iid)
#'   }
#'
#' }
#'
#' # outside of function ----------
#'
#' list_species_iid <- c ("alpaca",
#'                        "cat",
#'                        "chicken",
#'                        "cow",
#'                        "dog",
#'                        "duck",
#'                        "fly",
#'                        "guinea_pig",
#'                        "horse",
#'                        "human",
#'                        "mouse",
#'                        "pig",
#'                        "rabbit",
#'                        "rat",
#'                        "sheep",
#'                        "turkey",
#'                        "worm",
#'                        "yeast")
