
# Set up NetworkHub cache -----

## Create a "directory" for the NetworkHub cache, where you can later on store and get the PPI data -------


#' initialize_NetworkHub
#'
#' @param nh_cachedir
#'
#' @return
#' @export
#'
#' @examples

# nh_cachedir = default directory name where the cache will be stored
initialize_NetworkHub <- function(nh_cachedir = "NetworkHub") {

  # define the cache directory, with tools::R_user_dir you can create a path to the cache, which = cache, return the path for caching purposes
  cache_dir <- tools::R_user_dir(nh_cachedir, which = "cache")

  bfc_nh <- BiocFileCache::BiocFileCache(cache_dir)

  return(bfc_nh)

}





