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
#' url_sdb <- urlmaker_stringdb(species = "Homo sapiens", version = "12.0", type == "PPI")
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

# build_graph() --------------------

#' build_graph()
#'
#' @param graph_data ppi data
#' @param anno_df dataframe of ppi data from corresponding database containing annotation (Ensembl, GeneSymbol, Entrez, Uniprot)
#' @param g graph object created with build_graph()
#' @param idtype_anno Character, specifies which identifier type is used in the
#' `anno_df`. Defaults to "gene_symbol"
#' @param data_source database
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#' @param add_info_nodes TRUE adding annotation infos to nodes
#'
#' @return my_graph
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' db_stringdb_df <- get_networkdata_stringdb(species = "Mus musculus",
#'                                            version = "12.0",
#'                                            cache = TRUE,
#'                                            get_annotation = FALSE,
#'                                            add_annotation = FALSE
#'                                           )
#'
#' db_stringdb_anno_df <- get_annotation_stringdb(species = "Mus musculus",
#'                                                version = "12.0",
#'                                                cache = TRUE
#'                                                )
#'
#' db_stringdb_ppi_anno_df <- add_annotation_stringdb(ppi_stringdb = db_stringdb_df,
#'                                                    anno_df = db_stringdb_anno_df,
#'                                                    species = "Mus musculus")
#'
#'
#' db_stringdb_igraph_object <- build_graph(graph_data = db_stringdb_ppi_anno_df,
#'                                          anno_df = db_stringdb_anno_df,
#'                                          idtype_anno = "gene_symbol",
#'                                          data_source = "stringdb",
#'                                          output_format = "igraph",
#'                                          min_score_threshold = "0.35",
#'                                          add_info_nodes = TRUE
#'                                          )
#' }

build_graph <- function(graph_data,
                        anno_df,
                        idtype_anno,
                        data_source = c("stringdb", "hint", "funcoup", "iid", "irefindex", "mint",
                                        "genemania", "huri", "stringdb",
                                        "pathwaycommons", "reactome", "innatedb", "biogrid", "intact"),
                        output_format = "igraph",
                        min_score_threshold = NULL,
                        add_info_nodes = TRUE) {

  # Match the provided data_source to the available options
  data_source <- match.arg(data_source)

  # Build the function name dynamically
  function_name <- paste0("build_graph_", data_source)

  # Get the function by its name
  build_function <- get(function_name, envir = environment())

  # Call the specific function using the arguments
  my_graph <- build_function(graph_data = graph_data,
                             output_format = output_format,
                             min_score_threshold = min_score_threshold)

  if (add_info_nodes) {
    my_graph <- add_info_from_dataframe_to_graph(anno_df = anno_df, g = my_graph, idtype_anno = "gene_symbol")
  }

  return(my_graph)
}

# add_info_from_dataframe_to_graph() ----------


#' add_info_from_dataframe_to_graph()
#'
#' @param anno_df dataframe of ppi data from corresponding database containing annotation (Ensembl, GeneSymbol, Entrez, Uniprot)
#' @param g graph object created with build_graph()
#' @param idtype_anno Character, specifies which identifier type is used in the
#' `anno_df`. Defaults to "gene_symbol"
#'
#' @return An igraph object
#'
#' @importFrom igraph set_vertex_attr V
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' db_stringdb_df <- get_networkdata_stringdb(species = "Mus musculus",
#'                                            version = "12.0",
#'                                            cache = TRUE,
#'                                            get_annotation = FALSE,
#'                                            add_annotation = FALSE
#'                                           )
#'
#' db_stringdb_anno_df <- get_annotation_stringdb(species = "Mus musculus",
#'                                                version = "12.0",
#'                                                cache = TRUE
#'                                                )
#'
#' db_stringdb_ppi_anno_df <- add_annotation_stringdb(ppi_stringdb = db_stringdb_df,
#'                                                    anno_df = db_stringdb_anno_df,
#'                                                    species = "Mus musculus")
#'
#'
#' db_stringdb_igraph_object <- build_graph_stringdb(graph_data = db_stringdb_ppi_anno_df,
#'                                                   output_format = "igraph",
#'                                                   min_score_threshold = "0.35"
#'                                                   )
#'
#' db_stringdb_igraph_object_info_added <- add_info_from_dataframe_to_graph(anno_df = db_stringdb_anno_df,
#'                                                                          g = db_stringdb_igraph_object,
#'                                                                          idtype_anno = "gene_symbol")
#'                                                                      )
#'
#'
#' db_stringdb_igraph_object_info_added
#' }
#'
#' # check that it worked:
#' ## igraph::vertex_attr_names(graph = db_stringdb_igraph_object_info_added)
#' ## igraph::V(db_stringdb_igraph_object_info_added)$attr_gene_symbol
#' # this should have much less NAs, see above - possible TODO to re-check
#'
add_info_from_dataframe_to_graph <- function(anno_df,
                                             g,
                                             idtype_anno = "gene_symbol") {
  #loop through annotation column names
  for (col_name in c("entrez_id", "gene_symbol", "uniprot_id", "ensembl_id")) {
    # if column name is in dataframe
    if (col_name %in% colnames(anno_df)) {

      #define attributin column names for igraph object
      attr_name <- paste0("attr_", col_name)
      #match values from graph and annotation dataframe
      matched_values <- anno_df[match(igraph::V(g)$name, anno_df[[idtype_anno]]), col_name, drop = TRUE]

      #add values to graph
      g <- igraph::set_vertex_attr(g, name = attr_name, value = matched_values)
    }
  }

  return(g)
}

# map_SE_on_graph() -------------

#' map_SE_on_graph
#'
#' @param se summerized experiment object (e.g. DE analysis resluts)
#' @param value_to_map select the value of the results you want to compare
#' @param de_name name of experiment results/ what to compare (e.g. "ifng_vs_naive")
#' @param id_in_graph annotation type of se dataframe rows/ what should be mapped in the graph
#' @param g igraph object
#'
#' @return g
#' @export
#'
#' @import visNetwork
#' @import igraph
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' db_stringdb_df <- get_networkdata_stringdb(species = "Homo sapiens",
#'                                            version = "12.0",
#'                                            cache = TRUE,
#'                                            get_annotation = FALSE,
#'                                            add_annotation = FALSE
#'                                           )
#'
#' db_stringdb_anno_df <- get_annotation_stringdb(species = "Homo sapiens",
#'                                                version = "12.0",
#'                                                cache = TRUE
#'                                                )
#'
#' db_stringdb_ppi_anno_df <- add_annotation_stringdb(ppi_stringdb = db_stringdb_df,
#'                                                    anno_df = db_stringdb_anno_df,
#'                                                    species = "Homo sapiens")
#'
#'
#' db_stringdb_igraph_object <- build_graph(graph_data = db_stringdb_ppi_anno_df,
#'                                          anno_df = db_stringdb_anno_df,
#'                                          idtype_anno = "gene_symbol",
#'                                          data_source = "stringdb",
#'                                          output_format = "igraph",
#'                                          min_score_threshold = "0.35",
#'                                          add_info_nodes = TRUE
#'                                         )
#' g <- db_stringdb_igraph_object
#'
#' library("macrophage") #demo dataset
#' library("DESeq2") #package to perform differential expression analysis
#' library("org.Hs.eg.db") #Annotation package for human genes with info about gene identifier, chromosome location, Entrez-Gene, GO-terms, ...
#' library("AnnotationDbi") #Annotation package
#'
#' data(gse, package = "macrophage")
#' dds_macrophage <- DESeqDataSet(gse, design = ~ line + condition)
#' rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
#' keep <- rowSums(counts(dds_macrophage) >= 10) >= 6
#' dds_macrophage <- dds_macrophage[keep, ]
#'
#' dds_macrophage <- DESeq(dds_macrophage)
#'
#' res_macrophage_IFNg_vs_naive <- results(dds_macrophage,
#'                                        contrast = c("condition", "IFNg", "naive"),
#'                                        lfcThreshold = 1,
#'                                        alpha = 0.05)
#'
#' res_de <- res_macrophage_IFNg_vs_naive
#' de_df <- mosdef::deresult_to_df(res_de)
#'
#' top_ones <- mapIds(org.Hs.eg.db, de_df$id[1:100], column = "SYMBOL", keytype = "ENSEMBL")
#' top_ones_narm <- top_ones[!is.na(top_ones)]
#'
#' top_ones_in_graph <- intersect(top_ones_narm, V(g)$name)
#'
#'
#' #SE object
#' dds <- dds_macrophage
#' de_res <- res_de
#' de_name <- "ifngVSnaive"
#'
#' #function to combine dds dataframe (SE) with the results de dataframe (df)
#' combine_dds_deres <- function(dds, de_res, de_name) {
#'    matched_ids <- match(rownames(res_de), rownames(dds))
#'    rowData(dds)[[paste0(de_name, "_log2FoldChange")]] <- NA
#'    rowData(dds)[[paste0(de_name, "_pvalue")]] <- NA
#'    rowData(dds)[[paste0(de_name, "_padj")]] <- NA
#'    rowData(dds)[[paste0(de_name, "_log2FoldChange")]][matched_ids] <- res_de$log2FoldChange
#'    rowData(dds)[[paste0(de_name, "_pvalue")]][matched_ids] <- res_de$pvalue
#'    rowData(dds)[[paste0(de_name, "_padj")]][matched_ids] <- res_de$padj
#'  return(dds)
#'  }
#' se_macrophage <- combine_dds_deres(dds = dds_macrophage, de_res = res_de, de_name = "ifng_vs_naive")
#'
#' # what is this for? #TODO
#' de_res_scrambles <- de_res
#' de_res_scrambles$log2FoldChange <- rnorm(17806, sd = 2)
#'
#' se_macrophage
#'
#' graph_se_macrophage_string <- map_SE_on_graph(se = se_macrophage,
#'                                               de_name = "ifng_vs_naive",
#'                                               value_to_map = "log2FoldChange",
#'                                               id_in_graph = "attr_ensembl_id",
#'                                               g = db_stringdb_igraph_object
#'                                               )
#'
#' #Extra: Creation of a visNetwork :)
#'
#' subg <- induced_subgraph(
#'    graph_se_macrophage_string,
#'    top_ones_in_graph
#' )
#'
#' visNetwork::visIgraph(subg) |>
#'   visNetwork::visEdges(color = list(color = "#88888888"))|>
#'     visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
#'                                                    degree = 1,
#'                                                    hover = TRUE),
#'                            nodesIdSelection = TRUE)
#'
#' }

map_SE_on_graph <- function(se,
                            de_name = "ifng_vs_naive",
                            value_to_map = c("log2FoldChange", "pvalue"),
                            id_in_graph = "attr_ensembl_id",
                            g
                            ) {

    # Initialize a vector with `NA` for all nodes
    values <- rep(NA, igraph::vcount(g))

    #take the names of the nodes from the graph
    node_names <- vertex_attr(g, id_in_graph)

    #match
    matched_names <- match(node_names, rownames(se))

    # extract value (e.g. log2foldchange) and save in vector
    column_name <- paste0(de_name, "_", value_to_map)
    values[!is.na(matched_names)] <- rowData(se)[[column_name]][matched_names[!is.na(matched_names)]]

    # colors <- ifelse(is.na(values), "grey",
    #                  ifelse(values > 0, "red", "blue"))
    #
    # # add the color as an attribute
    g <- igraph::set_vertex_attr(g, value_to_map, value = values)


    mypal <- rev(scales::alpha(
      colorRampPalette(RColorBrewer::brewer.pal(name = "RdBu", 11))(50), 0.4
    ))

    max_de <- max(abs(range(
      vertex_attr(g, value_to_map), na.rm = TRUE
    )))

    V(g)$color <- mosdef::map_to_color(vertex_attr(g, value_to_map), mypal, limits = c(-max_de, max_de))

    rank_gs <- rank(V(g)$name)
    g <- permute(g, rank_gs)

    # title for tooltips
    g <- set_vertex_attr(g, name = "title", value = NA)


    # vertex_titles <- paste("Gene Symbol:", vertex_attr(g, "attr_gene_symbol"),
    #                         "Entrez ID:", vertex_attr(g, "attr_entrez_id"),
    #                         "Uniprot ID:", vertex_attr(g, "attr_uniprot_id"),
    #                         "log2FC:", vertex_attr(g, "log2FoldChange")
    #                       )
    #
    #
    # g <- set_vertex_attr(g, name = "title", value = vertex_titles)

    vertex_title_genesymbol <- paste("Gene Symbol:", vertex_attr(g, "attr_gene_symbol"))
    vertex_title_entrez     <- paste("Entrez ID:", vertex_attr(g, "attr_entrez_id"))
    vertex_title_uniprot    <- paste("Uniprot ID:", vertex_attr(g, "attr_uniprot_id"))
    vertex_title_log2fc <- paste("log2FC:", round(vertex_attr(g, "log2FoldChange"), 4))

    vertex_titles <- paste0(vertex_title_genesymbol, "<br>",
                           vertex_title_entrez, "<br>",
                           vertex_title_uniprot, "<br>",
                           vertex_title_log2fc)

    #TODO title with line break
    g <- set_vertex_attr(g, name = "title", value = vertex_titles)

    rank_gs <- rank(V(g)$name)
    g <- permute.vertices(g, rank_gs)

    return(g)
}
