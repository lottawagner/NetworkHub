# Whats inside?
# 1. a function to initialize_NetworkHub cache on the corresponding machine
# 2. a function to cache_NetworkHub files from dbs
# 3. a function to get_NetworkHub files out of the cache

# initialize_NetworkHub() -------
#' initialize_NetworkHub()
#'
#' @param nh_cachedir default directory name where the cache will be stored
#'
#' @return bfc_nh variabe assigned to the chache
#'
#' @importFrom BiocFileCache BiocFileCache
#'
#' @export
#'
#' @examples
#' library("BiocFileCache") # load the library of the package
#' bfc_nh <- initialize_NetworkHub() # use this function to create the cache called "bfc_nh"
#' bfccache(bfc_nh) # get the path to the cache
#' length(bfc_nh) # how many entries?
#' bfcinfo(bfc_nh) # what are the column names of the cache
#' # bfcremove(bfc_nh, c("BFC10", "BFC11", "BFC12", "BFC13"))
#' # nh_cachedir = default directory name where the cache will be stored
initialize_NetworkHub <- function(nh_cachedir = "NetworkHub") {

  # define the cache directory, with tools::R_user_dir you can create a path to the cache, which = cache, return the path for caching purposes
  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache")

  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir) # creats/uploads objects to cache

  return(bfc_nh)

}


# cache_NetworkHub() ----

#' cache_NetworkHub()
#'
#' @param rname resource name in the cache (local)
#' @param fpath filepath (external)
#' @param nh_cachedir default directory name where the cache will be stored
#' @param ... further arguments passed to or from other methods
#'
#' @return rpath returns the resource path of the cached file
#'
#' @importFrom BiocFileCache BiocFileCache bfcquery bfccount bfcadd
#'
#' @export
#'
#' @examples
#' # for example, retrieve something from stringDB
#'
#' url_sdb <- urlmaker_stringdb("PPI", "Homo sapiens", "12.0")
#' url_sdb
#'
#' cache_NetworkHub(
#'   rname = "STRINGDB_Homo sapiens_v12.0",
#'   fpath = url_sdb
#' )

cache_NetworkHub <- function(rname, # ressource name
                             fpath, # filepath
                             nh_cachedir = "NetworkHub", # name of cache
                             ...) { # additional arguments transfered to bfcadd
  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache") # tells the function to use the previously defined directory nh_cachedir
  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir) # creats/uploads objects to cache

  # check if filepath (fpath) is already being tracked

  nh_query <- BiocFileCache::bfcquery(bfc_nh, fpath)
  if (BiocFileCache::bfccount(nh_query) == 0) {
    rpath <- BiocFileCache::bfcadd(bfc_nh, rname, fpath, ...) # add it to the cache by identifing the ressource with rname & fpath and ... (additional arguments that are passed to cache)
  }
  if (BiocFileCache::bfccount(nh_query) == 1){ # == 1 -> fpath once in cache?
    rpath <- nh_query$rpath # use the already existing ressource path
  }
  if (BiocFileCache::bfccount(nh_query) > 1) {
    warning("WARNING - more than one version with the same fpath contained in bfc_nh", nh_query$fpath)
    rpath <- nh_query$rpath[1]
  }

  return(rpath) # returns the resource path of the cached file
}

# fetch_NetworkHub() ----

#' fetch_NetworkHub()
#'
#' @param rname resource name in the cache (local)
#' @param update default parameter set to TRUE (up to date version of the resource)
#' @param nh_cachedir default directory name where the cache will be stored
#' @param ... further arguments passed to or from other methods
#'
#' @return res_nh returns rpath to res_nh or NULL if no entry
#' @export
#'
#' @examples
#' fetch_NetworkHub(rname = "STRINGDB_Homo sapiens_v12.0")
#'
fetch_NetworkHub <- function(rname, # resourcename
                           update = TRUE, # up to date version of the resource
                           nh_cachedir = "NetworkHub", # name of cache
                           ...) { # additional arguments transferred to bfcadd

  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache") # tells the function to use the previously defined directory nh_cachedir
  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir) # creates/uploads objects to cache

  nh_query <- BiocFileCache::bfcquery(bfc_nh, rname, exact = TRUE) # search for entry rname and saves in variable nh_query

  # is there already a cached version?
  res_nh <- NULL # the result should be saved in res_nh

  if (BiocFileCache::bfccount(nh_query)) { # if there is already the file
    rid <- nh_query$rid # takes the ressource ID (rid) from the entry

    nu <- FALSE # is there a new update necessary?
    if (update) { #TRUE?
      nu <- BiocFileCache::bfcneedsupdate(bfc_nh, rid) # nu (newest update) is defined to the rid in the cache
    }
    if (!isFALSE(nu)) {
      BiocFileCache::bfcdownload(bfc_nh, rid, ask = FALSE, ...) #download of the newest version of the entry
    }

    message("Using cached version from ", nh_query$create_time) # tells you from which date this version is
    res_nh <- BiocFileCache::bfcrpath(bfc_nh, rname) # rpath of updated entry is saved in variable res_nh
  }
if (is.null(res_nh)) {
  message("No record found in NetworkHub!") # if res_nh is NULL (no file in cache) -> tells you nothing found
}
  return(res_nh) # returns rpath to res_nh or NULL if no entry
}

## build_graph() allgemeine Funktion? --------------------

#' #' Title
#' #'
#' #' @param graph_data data describing the nodes and edges of the network
#' #' @param data_source database (e.g.: "stringdb")
#' #' @param output_format vector that specifies possible output formats ("igraph", "grpahnel", "sif")
#' #' @param min_score_threshold optional threshold to filter edges based on their combination score
#' #' @param subset_nodes optional vector of nodes to extract only a part of the network
#' #'
#' #' @return my_graph object
#' #' @export
#' #'
#' #' @examples
#' #'
#' #' TODO
#'
#' build_graph <- function(graph_data, # data describing the nodes and edges of the network
#'                         data_source = c("stringdb", "hint", "funcoup", "iid", "irefindex", "mint", "genemania", "reactome", "cpdb", "huri", "matrixdb", "pathwaycommons", "innatedb", "biogrid", "corum", "intact"),
#'                         output_format = c("igraph", "grpahnel", "sif"), # Vector that specifies the desired output formats
#'                         min_score_threshold, # optional threshold for filtering the edges based on a combination score
#'                         subset_nodes){ # optional subset of nodes that can be used to create a subgraph
#'
#'   # define the data_source and the corresponding handling for the graph
#'   if (data_source == "stringdb") {
#'     my_graph <- build_graph_stringdb(
#'       graph_data = graph_data,
#'       output_format = output_format,
#'       min_score_threshold = min_score_threshold,
#'       subset_nodes = subset_nodes
#'     )
#'
#'     # still needs to be defined
#'   } else {
#'     #TODO once we worked on the other databases
#'   }
#'
#'   #return my_graph
#'   return (my_graph)
#' }

