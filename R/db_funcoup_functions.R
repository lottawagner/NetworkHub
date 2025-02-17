# get_networkdata_funcoup() -----------

#' get_networkdata_funcoup()
#'
#' @param species  from which species does the data come from c( "A.thaliana", "B.subtilis", "B.taurus", "C.elegans","C.familiaris", "C.intestinalis", "D.melanogatser", "D.rerio", "E.coli", "G.gallus", "H.sapiens", "M.jannaschii", "M.musculus", "O.sativa", "P.falciparum", "R.norvegicus", "S.cerevisiae", "S.pombe", "S.scrofa", "S.solfataricus")
#' @param version version of the data files in funcoup
#' @param type compact or full file
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return A data frame containing the PPI information in tabular format
#'
#' @importFrom vroom vroom
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_funcoup_df <- get_networkdata_funcoup(
#'   species = "H.sapiens",
#'   version = "6.0",
#'   type = "compact",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#' db_funcoup_df
#' }
#'
get_networkdata_funcoup <- function(species = "H.sapiens",
                                    version = "6.0",
                                    type = c("compact", "full"),
                                    cache = TRUE,
                                    get_annotation = TRUE,
                                    add_annotation = TRUE,
                                    ...) {
  # list species is actualized for version funcoup v.5.0 build 2020-09
  # UPDATEVERSION

  # check that the value for species is listed in funcoup

  if (!(species %in% list_species_funcoup)) { # if species is not in the list
    stop(
      "Species not found as specified by FunCoup,",
      "please check some valid entries by running `list_species_funcoup`"
    ) # stop function and print
  }

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "funcoup_",
    species,
    "_v",
    version,
    "_",
    type
  ) # definition of the resource name

  # TODO problem with caching what is happening? it is always downloading again (3x)
  if (cache) {
    # tries to fetch from the cache
    message("Fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("File not in cache, downloading to cache...")
    # buildup from funcoup url
    funcoup_url <-
      urlmaker_funcoup(
        species = species,
        version = version,
        type = type
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_funcoup

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = funcoup_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppi_funcoup <- vroom::vroom(network_file)
  # ppi_funcoup <- head(read.delim(network_file, sep = "\t"))

  message(dim(ppi_funcoup))
  # colnames(ppi_funcoup)

  colnames(ppi_funcoup)[colnames(ppi_funcoup) == "0:ProteinA"] <- "Uniprot_A"
  colnames(ppi_funcoup)[colnames(ppi_funcoup) == "1:ProteinB"] <- "Uniprot_B"

  annotation_db <-
    funcoup_db_annotations$anno_db_funcoup[match(species, funcoup_db_annotations$species)]


  if (get_annotation && is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_funcoup(..., add_annotation = FALSE, ...)\n"
    )
    return(ppi_funcoup)
  }

  if (get_annotation && !is.na(annotation_db)) {
    db_funcoup_anno_df <- get_annotation_funcoup(
      ppi_funcoup = ppi_funcoup,
      species = species
    )

    message("...created annotation dataframe")

    if (add_annotation) {
      db_funcoup_ppi_anno_df <- add_annotation_funcoup(
        anno_df = db_funcoup_anno_df,
        ppi_funcoup = ppi_funcoup,
        species = species
      )

      message("...added annotation from *db_funcoup_anno_df* to *db_funcoup_ppi_anno_df*")

      return(db_funcoup_ppi_anno_df)
    }

    if (!add_annotation) {
      return(db_funcoup_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation) {
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppi_funcoup)
}


# outside of function ----------

list_species_funcoup <- c(
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

list_db_annotationdbi_funcoup <- c(
  "org.At.tair.db",
  NA,
  "org.Bt.eg.db",
  "org.Ce.eg.db",
  "org.Cf.eg.db",
  NA,
  "org.Dm.eg.db",
  "org.Dr.eg.db",
  "org.EcK12.eg.db",
  "org.Gg.eg.db",
  "org.Hs.eg.db",
  NA,
  "org.Mm.eg.db",
  NA,
  NA,
  "org.Rn.eg.db",
  "org.Sc.sgd.db",
  NA,
  "org.Ss.eg.db",
  NA
)


funcoup_db_annotations <- data.frame(
  species = list_species_funcoup,
  anno_db_funcoup = list_db_annotationdbi_funcoup,
  row.names = list_species_funcoup
)


# get_annotation_funcoup() --------

#' get_annotation_funcoup ()
#'
#' @param species from which species does the data come from c( "A.thaliana", "B.subtilis", "B.taurus", "C.elegans","C.familiaris", "C.intestinalis", "D.melanogatser", "D.rerio", "E.coli", "G.gallus", "H.sapiens", "M.jannaschii", "M.musculus", "O.sativa", "P.falciparum", "R.norvegicus", "S.cerevisiae", "S.pombe", "S.scrofa", "S.solfataricus")
#' @param ppi_funcoup variable defined by ppi_funcoup in get_networkdata_funcoup()
#'
#' @importFrom AnnotationDbi mapIds
#'
#' @import org.At.tair.db
#' @import org.Bt.eg.db
#' @import org.Ce.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.Dr.eg.db
#' @import org.EcK12.eg.db
#' @import org.Gg.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#' @import org.Ss.eg.db
#'
#' @return anno_df
#' @return ppi_funcoup
#'
#' @export
#'
#' @examples
#' \dontrun{
#' db_funcoup_df <- get_annotation_funcoup(ppi_funcoup,
#'   species = "H.sapiens"
#' )
#' }
get_annotation_funcoup <- function(ppi_funcoup,
                                   species) {
  # find database on corresponding species

  if (!(species %in% list_species_funcoup)) { # if species is not in the list
    stop(
      "Species not found as specified by FunCoup,",
      "please check some valid entries by running `list_species_funcoup`"
    ) # stop function and print
  }


  annotation_db <-
    funcoup_db_annotations$anno_db_funcoup[match(species, funcoup_db_annotations$species)]

  if (!is.na(annotation_db)) {
    all_prot_ids <- unique(c(ppi_funcoup$Uniprot_A, ppi_funcoup$Uniprot_B))

    anno_df <- data.frame(
      uniprot_id = all_prot_ids,
      gene_symbol = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "SYMBOL", multiVals = "first"
      ),
      ensembl_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENSEMBL", multiVals = "first"
      ),
      entrez_id = mapIds(
        get(annotation_db),
        keys = all_prot_ids, keytype = "UNIPROT", column = "ENTREZID", multiVals = "first"
      ),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    message(
      "Annotation database for the species is not implemented yet.\n",
      "Next time define add_annotation in get_networkdata_funcoup(..., add_annotation = FALSE, ...)\n"
    )
    return(ppi_funcoup)
  }
}

# add_annotation_funcoup () -------------
#' add_annotation_funcoup ()
#'
#' @param anno_df annotation dataframe (for corresponding species in funcoup)
#' @param ppi_funcoup variable defined by ppi_funcoup in get_networkdata_funcoup()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_funcoup with annotation columns for each interactor (for corresponding species in funcoup)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_funcoup_df <- get_networkdata_funcoup(
#'   species = "H.sapiens",
#'   version = "5.0",
#'   cache = TRUE,
#'   get_annotation = FALSE,
#'   add_annotation = FALSE
#' )
#'
#' db_funcoup_anno_df <- get_annotation_funcoup(
#'   ppi_funcoup = db_funcoup_df,
#'   species = "H.sapiens"
#' )
#'
#' db_funcoup_ppi_anno_df <- add_annotation_funcoup(
#'   ppi_funcoup = db_funcoup_df,
#'   anno_df = db_funcoup_anno_df,
#'   species = "H.sapiens"
#' )
#' }
add_annotation_funcoup <- function(ppi_funcoup,
                                   anno_df,
                                   species) {
  ppi_funcoup <- ppi_funcoup

  # adding GeneSymbol
  ppi_funcoup$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppi_funcoup$Uniprot_A, anno_df$uniprot_id)]
  ppi_funcoup$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppi_funcoup$Uniprot_B, anno_df$uniprot_id)]
  # adding Uniprot
  ppi_funcoup$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_funcoup$Uniprot_A, anno_df$uniprot_id)]
  ppi_funcoup$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_funcoup$Uniprot_B, anno_df$uniprot_id)]
  # adding Entrez
  ppi_funcoup$Entrez_A <-
    anno_df$entrez_id[match(ppi_funcoup$Uniprot_A, anno_df$uniprot_id)]
  ppi_funcoup$Entrez_B <-
    anno_df$entrez_id[match(ppi_funcoup$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_funcoup)
}



# build_graph_funcoup() -----

#' build_graph_funcoup()
#'
#' @param graph_data ppi data from funcoup
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#'
#' @importFrom igraph graph.data.frame simplify
#'
#' @return my_graph
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_funcoup_df <- get_networkdata_funcoup(
#'   species = "H.sapiens",
#'   version = "6.0",
#'   type = "compact"
#' )
#'
#' db_funcoup_graph <- build_graph_funcoup(
#'   graph_data = db_funcoup_df,
#'   output_format = "igraph",
#'   min_score_threshold = "0.2"
#' )
#' db_funcoup_graph # list of 17004
#' }
#'
build_graph_funcoup <- function(graph_data,
                                output_format = "igraph",
                                min_score_threshold = NULL) {
  # check on the clumns in your ppi data file
  colnames(graph_data)
  hist(graph_data$`#0:PFC`, breaks = 50)

  # select ppi data >= minimal score
  if (!is.null(min_score_threshold)) {
    graph_data_processed <- graph_data[graph_data$`#0:PFC` >= min_score_threshold, ]
  } else {
    graph_data_processed <- graph_data
  }

  # check on dimension (amount of rows)
  dim(graph_data)
  dim(graph_data_processed)

  edges <- data.frame(
    from = graph_data_processed$GeneSymbol_A,
    to = graph_data_processed$GeneSymbol_B
  )

  # Create unique nodes (combine both GeneSymbol columns)
  nodes <- data.frame(
    id = unique(c(
      graph_data_processed$GeneSymbol_A,
      graph_data_processed$GeneSymbol_B
    )),
    label = unique(c(
      graph_data_processed$GeneSymbol_A,
      graph_data_processed$GeneSymbol_B
    ))
  )

  # If output format is igraph, return the igraph object
  if (output_format == "igraph") {
    whole_graph <- igraph::graph.data.frame(d = edges, directed = FALSE)
    my_graph <- igraph::simplify(whole_graph)
    return(my_graph)
  }
  # simplify by avoiding multiple entries?
  ## could make it smaller and easier to handle, without losing too much/at all in info
}
