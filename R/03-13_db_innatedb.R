# get_networkdata_innatedb() -----------

#' get_networkdata_innatedb()
#'
#' @param species  default value = "taxid:9606(Homo sapiens)" - from which species does the data come from
#' @param version version of the data files in innatedb
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param add_annotation expanding the dataframe with six columns (GeneSymbol,Entrez_ID and Ensembl_ID)
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_innatedb
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_innatedb_df <- get_networkdata_innatedb(species = "taxid:9606(Human)",
#'                                            version = "5.4"
#'                                            )
#' db_innatedb_df
#' }

get_networkdata_innatedb <- function(species,
                                     version = "5.4",
                                     cache = TRUE,
                                     add_annotation = TRUE,
                                     ...) {

  # list species is actualized for version innatedb "current" (sept 2024)
  # UPDATEVERSION

  # check that the value for species is listed in innatedb

  if (!(species %in% list_species_innatedb)) { # if species is not in the list
    stop("Species not found as specified by innatedb,",
         "please check some valid entries by running `list_species_innatedb`") # stop function and print
  }


  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0("innatedb_v_",
                  version
                  )

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" innatedb url
    innatedb_url <-
      urlmaker_innatedb(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_innatedb

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = innatedb_url
    )

  }

  # read in the resource, whether cached or freshly downloaded
  ppis_innatedb <- vroom::vroom(network_file)
  #ppis_innatedb <- head(read.delim(network_file, sep = " "))

  #If you want to check all species contained in the file
  #species_A <- unique(ppis_innatedb$ncbi_taxid_A)
  #species_B <- unique(ppis_innatedbk$ncbi_taxid_B)
  #info_species_innatedb <- union(species_A, species_B)

  #Filter for species

  ppis_innatedb_filtered <- ppis_innatedb[(ppis_innatedb$ncbi_taxid_A == species ) &
                                          (ppis_innatedb$ncbi_taxid_B == species ),
                                          ]

  # rename columns

  #Uniprot
  colnames(ppis_innatedb_filtered)[colnames(ppis_innatedb_filtered) == "alt_identifier_A"] <- "Ensembl_A"
  colnames(ppis_innatedb_filtered)[colnames(ppis_innatedb_filtered) == "alt_identifier_B"] <- "Ensembl_B"

  #extract Ensembl_id
  ppis_innatedb_filtered$Ensembl_A <- str_extract(ppis_innatedb_filtered$Ensembl_A, "ensembl:([A-Z0-9]+)")
  ppis_innatedb_filtered$Ensembl_A <- gsub("ensembl:", "", ppis_innatedb_filtered$Ensembl_A)

  ppis_innatedb_filtered$Ensembl_B <- str_extract(ppis_innatedb_filtered$Ensembl_B, "ensembl:([A-Z0-9]+)")
  ppis_innatedb_filtered$Ensembl_B <- gsub("ensembl:", "", ppis_innatedb_filtered$Ensembl_B)


  #Confidence score
  ppis_innatedb_filtered$confidence_score <- str_extract(ppis_innatedb_filtered$confidence_score, "lpr:([0-9\\.]+)")
  ppis_innatedb_filtered$confidence_score <- gsub("lpr:", "", ppis_innatedb_filtered$confidence_score)


  if (add_annotation) {

    ppi_innatedb_filtered_annotated <- annotation_innatedb(ppi_innatedb = ppis_innatedb_filtered,
                                                           species = species,
                                                           version = version)
    return(ppi_innatedb_filtered_annotated)
  }

  if (!add_annotation) {
    return(ppis_innatedb_filtered)
  }
}


# Outside function ----------


#SPECIES NetworkHub uses from innatedb
list_species_innatedb <- c("taxid:9606(Human)",
                           "taxid:10090(Mouse)"
                          )


list_db_annotationdbi_innatedb <- c("org.Hs.eg.db",
                                    "org.Mm.eg.db"
                                   )


innatedb_db_annotations <- data.frame(species = list_species_innatedb,
                                      anno_db_innatedb = list_db_annotationdbi_innatedb,
                                      row.names = list_species_innatedb
)


# annotation_innatedb() --------

#' annotation_innatedb ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in innatedb
#' @param ppi_innatedb variable defined by ppis_innatedb in get_networkdata_innatedb()
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db

#'
#'@return ppis_innatedb_annotated
#'
#'@export
#'
#'
#' @examples
#' #\dontrun{
#' #annotation_innatedb <- annotation_innatedb(ppi_innatedb,
#'                                             species = "taxid:9606(Homo sapiens)",
#'                                             version = "current")
#' #annotation_innatedb
#' #}

annotation_innatedb <- function(ppi_innatedb,
                                species,
                                version)
{
  # find database on corresponding species

  if (!(species %in% list_species_innatedb)) { # if species is not in the list
    stop("Species not found as specified by innatedb,",
         "please check some valid entries by running `list_species_innatedb`") # stop function and print
  }

  annotation_db <-
    innatedb_db_annotations$anno_db_innatedb[match(species, innatedb_db_annotations$species)]


  all_prot_ids <- unique(na.omit(c(ppi_innatedb$Ensembl_A, ppi_innatedb$Ensembl_B)))

  anno_df <- data.frame(ensembl_id = all_prot_ids,
                        gene_symbol <- mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "ENSEMBL",
                                              column = "SYMBOL"
                        ),
                        uniprot_id <- mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "ENSEMBL",
                                             column = "UNIPROT"
                        ),
                        entrez_ids <- mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "ENSEMBL",
                                             column = "ENTREZID"
                        ),
                        row.names = all_prot_ids
  )

  ppis_innatedb_annotated <- ppi_innatedb

  #adding GeneSymbol
  ppis_innatedb_annotated$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppis_innatedb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_innatedb_annotated$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppis_innatedb_annotated$Ensembl_B, anno_df$ensembl_id)]

  #adding Uniprot
  ppis_innatedb_annotated$Uniprot_A <-
    anno_df$uniprot_id[match(ppis_innatedb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_innatedb_annotated$Uniprot_B <-
    anno_df$uniprot_id[match(ppis_innatedb_annotated$Ensembl_B, anno_df$ensembl_id)]

  #adding Entrez
  ppis_innatedb_annotated$Entrez_A <-
    anno_df$entrez_id[match(ppis_innatedb_annotated$Ensembl_A, anno_df$ensembl_id)]
  ppis_innatedb_annotated$Entrez_B <-
    anno_df$entrez_id[match(ppis_innatedb_annotated$Ensembl_B, anno_df$ensembl_id)]

  return(ppis_innatedb_annotated)
}



# build_graph_innatedb() -----

#' build_graph_innatedb()
#'
#' @param graph_data ppi data from innatedb
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident": lpr score (lowest PMID re-use)
#'
#' @importFrom igraph graph.data.frame simplify
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_innatedb_df <- get_networkdata_innatedb(species = "taxid:9606(Human)",
#'                                            version = "5.4"
#'                                            )
#'
#' db_innatedb_graph <- build_graph_innatedb(graph_data = db_innatedb_df,
#'                                           output_format = "igraph",
#'                                           min_score_threshold = "20" )
#' db_innatedb_graph #list of 690 (score = 10) 1835 (score = 20)
#' }
#'
#'

build_graph_innatedb <- function (graph_data,
                                  output_format = "igraph",
                                  min_score_threshold = NULL){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$confidence_score <- as.numeric(graph_data$confidence_score)

  # Erstelle das Histogramm mit 50 bins (breaks)
  hist(graph_data$confidence_score, breaks = 50)

  #select ppi data >= minimal score
  if (!is.null(min_score_threshold)){
    graph_data_processed <- graph_data[graph_data$confidence_score <= min_score_threshold, ]
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
  # simplify by avoiding multiple entries?
  ## could make it smaller and easier to handle, without losing too much/at all in info

}








