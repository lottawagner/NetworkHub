

## get_networkdata_stringdb() ---------

#' get_networkdata_stringdb()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in stringdb
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_stringdb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_string_df <- get_networkdata_stringdb(species = "Mus musculus",
#'                                          version = "12.0"
#'                                         )
#' db_string_df
#' }
#'
get_networkdata_stringdb <- function(species,
                                     version,
                                     cache = TRUE,
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
  #Uniprot
  colnames(ppis_stringdb)[colnames(ppis_stringdb) == "protein1"] <- "Ensembl_A"
  colnames(ppis_stringdb)[colnames(ppis_stringdb) == "protein2"] <- "Ensembl_B"
  ppis_stringdb$Ensembl_A <- str_extract(ppis_stringdb$Ensembl_A, "(EN[S][A-Z0-9]+)")
  ppis_stringdb$Ensembl_B <- str_extract(ppis_stringdb$Ensembl_B, "(EN[S][A-Z0-9]+)")

  if (add_annotation) {

    if (!(species %in% list_common_species_stringdb)) { # if species is not in the list
      stop("Species not in `list_common_species_stringdb`!",
           "Annotation for this species is not provided") # stop function and print
    }

    ppis_stringdb_annotated <- annotation_stringdb(ppi_stringdb = ppis_stringdb,
                                                        species = species,
                                                        version = version)
    message("...added annotation :)")
    return(ppis_stringdb_annotated)
  }

  if (!add_annotation){
    return(ppis_stringdb)
  }

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


# annotation_stringdb() --------

#' annotation_stringdb ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in stringdb
#' @param ppi_stringdb variable defined by ppis_stringdb in get_networkdata_stringdb()
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
#' @import org.At.tair.db
#' @import org.Bt.eg.db
#' @import org.Ce.eg.db
#' @import org.Cf.eg.db
#' @import org.Dm.eg.db
#' @import org.EcK12.eg.db
#' @import org.Gg.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Rn.eg.db
#' @import org.Sc.sgd.db
#' @import org.Ss.eg.db
#' @import org.Xl.eg.db
#'
#'@return ppis_stringdb_annotated
#'
#'@export
#'
#'
#' @examples
#' #\dontrun{
#' #annotation_stringdb <- annotation_stringdb(ppi_stringdb, species = "taxid:9606(Homo sapiens)", version = "current")
#' #annotation_stringdb
#' #}

annotation_stringdb <- function(ppi_stringdb,
                              species,
                              version){

  annotation_db <-
    stringdb_db_annotations$anno_db_stringdb[match(species, stringdb_db_annotations$species)]

  #create a list that contains all uniprot ids (but not NA)
  all_prot_ids <- unique(c(ppi_stringdb$Ensembl_A, ppi_stringdb$Ensembl_B))
  all_prot_ids <- gsub("\\.[0-9]+$", "", all_prot_ids)

  anno_df <- data.frame(ensembl_id = all_prot_ids,
                        gene_symbol = mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "ENSEMBLPROT",
                                              column = "SYMBOL"
                        ),
                        uniprot_id = mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "ENSEMBLPROT",
                                             column = "UNIPROT"
                        ),
                        entrez_ids = mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "ENSEMBLPROT",
                                             column = "ENTREZID"
                        ),
                        row.names = all_prot_ids
  )

  ppis_stringdb_annotated <- ppi_stringdb

  #adding GeneSymbol
  ppis_stringdb_annotated$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppis_stringdb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_stringdb_annotated$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppis_stringdb_annotated$Ensembl_B, anno_df$ensembl_id)]

  #adding Uniprot
  ppis_stringdb_annotated$Ensembl_A <-
    anno_df$ensembl_id[match(ppis_stringdb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_stringdb_annotated$Ensembl_B <-
    anno_df$ensembl_id[match(ppis_stringdb_annotated$Ensembl_B, anno_df$ensembl_id)]

  #adding Entrez
  ppis_stringdb_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_stringdb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_stringdb_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_stringdb_annotated$Ensembl_B, anno_df$ensembl_id)]

  return(ppis_stringdb_annotated)
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
#'                                        version =  "12.0")
#'
#' db_stringdb_graph <- build_graph_stringdb(graph_data = db_stringdb_df,
#'                                         output_format = "igraph",
#'                                         min_score_threshold = "0.35")
#' db_stringdb_graph #list of 17490
#' }
#'
#'


build_graph_stringdb <- function (graph_data,
                                output_format = "igraph",
                                min_score_threshold = NULL ){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$combined_score <- as.numeric(graph_data$combined_score)

  # Erstelle das Histogramm mit 50 bins (breaks)
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
