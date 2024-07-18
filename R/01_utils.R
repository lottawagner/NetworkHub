
# For what do we need the initialize step ?

initialize_NetworkHub <- function(nh_cachedir = "NetworkHub") {
  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache")
  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir)

  return(bfc_nh)

}

bfccache(bfc_nh)

#' @examples
#' bfc_nh <- initialize_NetworkHub ()
#' bfccache(bfc_nh)
#' length(bfc_nh)
#' show(bfc_nh)
#' bfcinfo(bfc_nh)
#'



# Create a "directory" for the NetworkHub cache, where you can later on store and get the PPI data -------

#' # Make sure to download the "BiocMananger" package before you start
#'
#' if (!"BiocManager" %in% rownames(installed.packages()))
#' install.packages("BiocManager")
#' BiocManager::install("BiocFileCache", dependencies=TRUE)
#'
#' library("BiocFileCache")


cache_NetworkHub <- function(rname, fpath, nh_cachedir = "NetworHub", ...) {
  cache.dir <- tools::R_user_dir(nh_cachedir, which = "cache")
  bfc_nh <- BiocFileCache::BiocFileCache(cache.dir)

  # check if fpath is being tracked
  nh_query <-  BiocFileCache::bfcquery(bfc_nh, fpath)
  if(BiocFileCache::bfccount(nh_query) == 0)
  {
    rpath <- BiocFileCache::bfcadd(bfc_nh, rname, fpath, ...)
  }
  else rpath <- nh_query$rpath
  return(rpath)
}


# get data out of the cache

get_NetworkHub <- function() {




}
