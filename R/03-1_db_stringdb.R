
# get_networkdata_stringdb() ---------

#' get_networkdata_stringdb()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in stringdb
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe , default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_stringdb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_string_df <- get_networkdata_stringdb(species = "Homo sapiens",
#'                                          version = "12.0",
#'                                          cache = TRUE,
#'                                          get_annotation = FALSE,
#'                                          add_annotation = FALSE
#'                                         )
#' db_string_df
#' }
#'
get_networkdata_stringdb <- function(species,
                                     version,
                                     cache = TRUE,
                                     get_annotation = TRUE,
                                     add_annotation = TRUE,
                                     ...){

  # matching for the species...
  species_id <- NULL # looking for the corresponding id in the following

  url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt",
                                  version)

  # read.delim the data from the species text file (columns separated using a delimiter)
  info_species <- read.delim(url(url_species_stringdb))
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
    match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]

  if (is.null(species_id))
    stop("No match found!")

  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "stringdb_",
    species,
    "_v",
    version,
    "_PPI"
  ) # definition of the  resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" stringdb url
    stringdb_url <-
      urlmaker_stringdb(
        type = "PPI",
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_stringdb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = stringdb_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_stringdb <- vroom::vroom(network_file)
  ## ppis_stringdb <- read.delim(network_file, sep = " ")


  # rename columns

  #Ensembl
  colnames(ppis_stringdb)[colnames(ppis_stringdb) == "protein1"] <- "Ensembl_Prot_A"
  colnames(ppis_stringdb)[colnames(ppis_stringdb) == "protein2"] <- "Ensembl_Prot_B"
  ppis_stringdb$Ensembl_Prot_A <- str_extract(ppis_stringdb$Ensembl_Prot_A, "(EN[S][A-Z0-9]+)")
  ppis_stringdb$Ensembl_Prot_B <- str_extract(ppis_stringdb$Ensembl_Prot_B, "(EN[S][A-Z0-9]+)")

  if (get_annotation) {

    if (!(species %in% list_common_species_stringdb)) { # if species is not in the list
      stop("Species not in `list_common_species_stringdb`!",
           "Annotation for this species is not provided") # stop function and print
    }

    db_stringdb_anno_df <- get_annotation_stringdb(species = species,
                                                   version = "12.0",
                                                   cache = TRUE
                                                  )

    message("...created annotation dataframe")


    if (add_annotation) {

      db_stringdb_ppi_anno_df <- add_annotation_stringdb(anno_df = db_stringdb_anno_df,
                                                         ppi_stringdb = ppis_stringdb,
                                                         species = species
      )
      message("...added annotation from *db_stringdb_anno_df* to *db_stringdb_ppi_anno_df*")

      return(db_stringdb_ppi_anno_df)

    }

    if (!add_annotation){
      return(db_stringdb_anno_df)
    }

  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }

  return(ppis_stringdb)

}



list_common_species_stringdb<- c("Arabidopsis thaliana",
                                "Bos taurus",
                                "Caenorhabditis elegans",
                                "Canis familiaris",
                                "Drosophila melanogaster",
                                "Escherichia coli",
                                "Gallus gallus",
                                "Homo sapiens",
                                "Mus musculus",
                                "Rattus norvegicus",
                                "Saccharomyces cerevisiae",
                                "Sus scrofa",
                                "Xenopus laevis")


list_db_annotationdbi_stringdb <- c("org.At.tair.db",
                                  "org.Bt.eg.db",
                                  "org.Ce.eg.db",
                                  "org.Cf.eg.db",
                                  "org.Dm.eg.db",
                                  "org.EcK12.eg.db",
                                  "org.Gg.eg.db",
                                  "org.Hs.eg.db",
                                  "org.Mm.eg.db",
                                  "org.Rn.eg.db",
                                  "org.Sc.sgd.db",
                                  "org.Ss.eg.db",
                                  "org.Xl.eg.db"
)


stringdb_db_annotations <- data.frame(species = list_common_species_stringdb,
                                    anno_db_stringdb = list_db_annotationdbi_stringdb,
                                    row.names = list_common_species_stringdb
)

# get_annotation_stringdb() ------

#' get_annotation_stringdb()
#'
#' @param species from which species does the data come from
#' @param version version of data files in STRING, default value = "12.0"
#' @param cache Default value set to TRUE (automatically checks if the data file is already stored in the cache)
#'
#' @return anno_df (for corresponding species in stringdb)
#' @export
#'
#' @examples
#'
#' db_stringdb_anno_df <- get_annotation_stringdb(species = "Homo sapiens",
#'                                                version = "12.0",
#'                                                cache = TRUE
#'                                                )
#'
#' # returns a dataframe containing 19699 rows and 5 columns
#'
#'
get_annotation_stringdb <- function(species,
                                    version,
                                    cache = TRUE) {
  # matching for the species...
  url_species_stringdb <- sprintf("https://stringdb-downloads.org/download/species.v%s.txt",
                                  version)

  # read.delim the data from the species text file (columns separated using a delimiter)
  info_species <- read.delim(url(url_species_stringdb))
  species_id <- info_species$X.taxon_id[ # in the column X.taxon_id we will find the taxon_id and assign it to the variable species_id
    match(species, info_species$official_name_NCBI) # matching the species to the corresponding entry in the info_species file column official_name_NCBI
  ]

  rname <- paste0(
    "STRINGDB_protein_info_",
    species,
    "_v",
    version
  )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    proteininfo_file <- fetch_NetworkHub(rname)
  }

  if (!cache | is.null(proteininfo_file)) {
    # retrieves the file for the first time
    message("Retrieving to cache...")
    # buildup from "base" STRINGDB url
    url_stringdb_info <-
      urlmaker_stringdb(
        type = "protein_info",
        species = species,
        version = version
      )

    proteininfo_file <- cache_NetworkHub(
      rname = rname,
      fpath = url_stringdb_info
    )
  }

  # read in the resource, whether cached or freshly downloaded
  stringdb_protein_info <- vroom::vroom(proteininfo_file)

  anno_df <- data.frame(
    protein_id = unique(stringdb_protein_info$`#string_protein_id`),
    row.names = unique(stringdb_protein_info$`#string_protein_id`)
  )

  # prefill with NAs
  anno_df$ensembl_id <- NA
  anno_df$gene_symbol <- NA
  anno_df$entrez_id <- NA
  anno_df$uniprot_id <- NA

  df_ensembl <- stringdb_protein_info[stringdb_protein_info$source == "Ensembl_gene", ]
  df_genesymbol <- stringdb_protein_info[stringdb_protein_info$source == "Ensembl_EntrezGene", ]
  df_entrez <- stringdb_protein_info[stringdb_protein_info$source == "Ensembl_HGNC_entrez_id", ]
  df_uniprot <- stringdb_protein_info[stringdb_protein_info$source == "Ensembl_UniProt", ]

  anno_df$ensembl_id <-
    df_ensembl$alias[match(anno_df$protein_id, df_ensembl$`#string_protein_id`)]
  anno_df$gene_symbol <-
    df_genesymbol$alias[match(anno_df$protein_id, df_genesymbol$`#string_protein_id`)]
  anno_df$entrez_id  <-
    df_entrez$alias[match(anno_df$protein_id, df_entrez$`#string_protein_id`)]
  anno_df$uniprot_id  <-
    df_uniprot$alias[match(anno_df$protein_id, df_uniprot$`#string_protein_id`)]

  anno_df$protein_id <- str_extract(anno_df$protein_id, "(EN[S][A-Z0-9]+)")

  return(anno_df)
}


# add_annotation_stringdb() --------

#' add_annotation_stringdb ()
#'
#' @param anno_df annotation dataframe (for corresponding species in stringdb)
#' @param ppi_stringdb variable defined by ppis_stringdb in get_networkdata_stringdb()
#' @param species  from which species does the data come from
#'
#'
#'@return ppi_stringdb with annotation columns for each interactor (for corresponding species in stringdb)
#'
#'@export
#'
#'
#' @examples
#' \dontrun{
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
#'                                                    species = "Homo sapiens"
#'                                                    )
#'
#'}

add_annotation_stringdb <- function(ppi_stringdb,
                                    anno_df,
                                    species
                                    ) {

  ppi_stringdb <- ppi_stringdb

  #adding Ensembl
  ppi_stringdb$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_stringdb$Ensembl_Prot_A, anno_df$protein_id)]
  ppi_stringdb$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_stringdb$Ensembl_Prot_B, anno_df$protein_id)]

  #adding GeneSymbol
  ppi_stringdb$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppi_stringdb$Ensembl_Prot_A, anno_df$protein_id)]
  ppi_stringdb$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppi_stringdb$Ensembl_Prot_B, anno_df$protein_id)]

  #adding Entrez
  ppi_stringdb$Entrez_A <-
    anno_df$entrez_id[match(ppi_stringdb$Ensembl_Prot_A, anno_df$protein_id)]
  ppi_stringdb$Entrez_B <-
    anno_df$entrez_id[match(ppi_stringdb$Ensembl_Prot_B, anno_df$protein_id)]

  #adding Uniprot
  ppi_stringdb$Uniprot_A <-
    anno_df$uniprot_id[match(ppi_stringdb$Ensembl_Prot_A, anno_df$protein_id)]
  ppi_stringdb$Uniprot_B <-
    anno_df$uniprot_id[match(ppi_stringdb$Ensembl_Prot_B, anno_df$protein_id)]

  return(ppi_stringdb)

}




# build_graph_stringdb() -----

#' build_graph_stringdb()
#'
#' @param graph_data ppi data from stringdb
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#'
#' @importFrom igraph graph.data.frame simplify
#' @importFrom graphics hist
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_stringdb_df <- get_networkdata_stringdb(species = "Homo sapiens",
#'                                            version = "12.0",
#'                                            cache = TRUE,
#'                                            get_annotation = TRUE,
#'                                            add_annotation = TRUE
#'                                           )
#'
#' db_stringdb_igraph_object <- build_graph_stringdb(graph_data = db_stringdb_df,
#'                                                   output_format = "igraph",
#'                                                   min_score_threshold = "0.35")
#' db_stringdb_graph #list of 17490
#' }


build_graph_stringdb <- function (graph_data,
                                  output_format = "igraph",
                                  min_score_threshold = NULL ){

  #check on the columns in your ppi data file
  colnames(graph_data)

  graph_data$combined_score <- as.numeric(graph_data$combined_score)

  # create hist with 50 breaks
  hist(graph_data$combined_score, breaks = 50)

  #select ppi data >= minimal score
  if (!is.null(min_score_threshold)){
    graph_data_processed <- graph_data[graph_data$combined_score >= min_score_threshold, ]
  } else {
    graph_data_processed <- graph_data
  }

  #check on dimension (amount of rows)
  dim(graph_data)
  dim(graph_data_processed)

  edges <- data.frame(from = graph_data_processed$GeneSymbol_A,
                      to = graph_data_processed$GeneSymbol_B)

  # Create unique nodes (combine both GeneSymbol columns)
  nodes <- data.frame(id = unique(c(graph_data_processed$GeneSymbol_A,
                                    graph_data_processed$GeneSymbol_B)),
                      label = unique(c(graph_data_processed$GeneSymbol_A,
                                       graph_data_processed$GeneSymbol_B)))

  # If output format is igraph, return the igraph object
  if (output_format == "igraph") {
    whole_graph <- igraph::graph.data.frame(d = edges, directed = FALSE)
    my_graph <- igraph::simplify(whole_graph)
    return(my_graph)
  }
}
