# Whats inside?
# 1. a function to initialize_NetworkHub cache on the corresponding machine
# 2. a function to cache_NetworkHub files from dbs
# 3. a function to get_NetworkHub files out of the cache

# Set up NetworkHub cache -----

## Create a "directory" for the NetworkHub cache, where you can later on store and get the PPI data -------

#QUESTION: Why does it make sense to create a function doing that? Is it easier for the user to just use the function?

#' initialize_NetworkHub function
#'
#' @param nh_cachedir
#'
#' @return
#' @export
#'
#' @examples
#' library("BiocFileCache") # load the library of the package
#' bfc_nh <- initialize_NetworkHub() # use this function to create the cache called "bfc_nh"
#' bfccache(bfc_nh) # get the path to the cache
#' length(bfc_nh) # how many entries?
#' bfcinfo(bfc_nh) # what are the column names of the cache

# nh_cachedir = default directory name where the cache will be stored

initialize_NetworkHub <- function(nh_cachedir = "NetworkHub") {

  # define the cache directory, with tools::R_user_dir you can create a path to the cache, which = cache, return the path for caching purposes
  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache")

  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir) # creats/uploads objects to cache

  return(bfc_nh)

}


# retrieve datafiles and store them in the cache ----

#' Title
#'
#' @param rname
#' @param fpath
#' @param nh_cachdir
#' @param ...
#'
#' @return
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

  nh_query <- BiocFileCache::bfcquery(bfc_nh, fpath) # is the fpath already in the cache?
  if (BiocFileCache::bfccount(nh_query) == 0) { # fpath is not in the cache already ?
    rpath <- BiocFileCache::bfcadd(bfc_nh, rname, fpath, ...) # add it to the cache by identifing the ressource with rname & fpath and ... (additional arguments that are passed to cache)
  } else { # ≠ 0 -> fpath already there?
    rpath <- nh_query$rpath # use the already existing ressource path
  }
  return(rpath) # retourns the ressource path of the cached file
}

# Get data from the cache ----

#' Title
#'
#' @param rname
#' @param update
#' @param nh_cachedir
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' get_NetworkHub(rname = "STRINGDB_Homo sapiens_v12.0")
#'
get_NetworkHub <- function(rname, # ressourcename
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


