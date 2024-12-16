# get_networkdata_hippie() ---------

#' get_networkdata_hippie()
#'
#' @param species default value = "Homo_sapiens", because this database only provides human data
#' @param version default value = "current", version of the database ... #UPDATEVERSION
#' @param cache default value = TRUE, (automatically checks if the data file is already stored in the cache)
#' @param get_annotation default value = TRUE, creation of an annotation dataframe , default value set to TRUE
#' @param add_annotation default value = TRUE, adding annotation to ppi datatframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_hippie
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_hippie_df <- get_networkdata_hippie(species = "Homo_sapiens",
#'                                        version = "current",
#'                                        cache = TRUE,
#'                                        get_annotation = FALSE,
#'                                        add_annotation = FALSE
#'                                        )
#' db_hippie_df
#' }
#'
get_networkdata_hippie <- function(species,
                                   version,
                                   cache = TRUE,
                                   get_annotation = TRUE,
                                   add_annotation = TRUE,
                                   ...){



  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "hippie_",
    species,
    "_v",
    version
  ) # definition of the  resource name

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" hippie url
    hippie_url <-
      urlmaker_hippie(
        species = species,
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_hippie

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = hippie_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_hippie <- vroom::vroom(network_file)
  ## ppis_hippie <- read.delim(network_file, sep = " ")


  # rename columns

  #Entrez
  colnames(ppis_hippie)[colnames(ppis_hippie) == "ID Interactor A"] <- "Entrez_A"
  ppis_hippie$Entrez_A <- gsub("entrez gene:", "", ppis_hippie$Entrez_A)
  colnames(ppis_hippie)[colnames(ppis_hippie) == "ID Interactor B"] <- "Entrez_B"
  ppis_hippie$Entrez_B <- gsub("entrez gene:", "", ppis_hippie$Entrez_B)

  #GeneSymbol
  colnames(ppis_hippie)[colnames(ppis_hippie) == "Gene Name Interactor A"] <- "GeneSymbol_A"
  colnames(ppis_hippie)[colnames(ppis_hippie) == "Gene Name Interactor B"] <- "GeneSymbol_B"

  if (get_annotation) {

    if (!(species %in% list_species_hippie)) { # if species is not in the list
      stop("Species not in `list_common_species_hippie`!",
           "Annotation for this species is not provided") # stop function and print
    }

    db_hippie_anno_df <- get_annotation_hippie(ppi_hippie = ppis_hippie,
                                               species = species,
                                               version = version
    )

    message("...created annotation dataframe")


    if (add_annotation) {

      db_hippie_ppi_anno_df <- add_annotation_hippie(anno_df = db_hippie_anno_df,
                                                     ppi_hippie = ppis_hippie,
                                                     species = species
      )
      message("...added annotation from *db_hippie_anno_df* to *db_hippie_ppi_anno_df*")

      return(db_hippie_ppi_anno_df)

    }

    if (!add_annotation){
      return(db_hippie_anno_df)
    }

  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }

  return(ppis_hippie)

}

# outside of functions --------


list_species_hippie<- c("Homo_sapiens")

list_db_annotationdbi_hippie <- c("org.Hs.eg.db")

hippie_db_annotations <- data.frame(species = list_species_hippie,
                                    anno_db_hippie = list_db_annotationdbi_hippie,
                                    row.names = list_species_hippie
                                    )

# get_annotation_hippie() ------

#TODO

#' get_annotation_hippie()
#'
#' @param ppi_hippie variable defined by ppis_hippie in get_networkdata_hippie()
#' @param species default value = "Homo_sapiens", from which species does the data come from
#' @param version default value = "current", version of data files in HIPPIE
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#'
#' @return anno_df (for corresponding species in hippie)
#' @export
#'
#' @examples
#'
#' db_hippie_df <- get_networkdata_hippie(species = "Homo_sapiens",
#'                                        version = "current",
#'                                        cache = TRUE,
#'                                        get_annotation = FALSE,
#'                                        add_annotation = FALSE
#'                                        )
#'
#' db_hippie_anno_df <- get_annotation_hippie(ppi_hippie = db_hippie_df,
#'                                            species = "Homo_sapiens",
#'                                            version = "current"
#'                                            )
#'
#' # returns a dataframe containing #TODO rows and #TODO columns
#'
#'
get_annotation_hippie <- function(ppi_hippie,
                                  species,
                                  version
                                  ) {
  annotation_db <-
    hippie_db_annotations$anno_db_hippie[match(species, hippie_db_annotations$species)]

  # unique(c(ppi_hippie$GeneSymbol_A, ppi_hippie$GeneSymbol_B))
  # damn, there is an NA entry inside the entries of hippie itself?!
  # table(is.na(unique(c(ppi_hippie$GeneSymbol_A, ppi_hippie$GeneSymbol_B))))

  if (!is.na(annotation_db)) {
    all_prot_ids <- na.omit(unique(c(ppi_hippie$GeneSymbol_A, ppi_hippie$GeneSymbol_B)))
    anno_df <- data.frame(
      genesymbol = all_prot_ids,
      ensembl_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "ENSEMBL"),
      uniprot_id = mapIds(
        get(annotation_db), keys = all_prot_ids, keytype = "SYMBOL", column = "UNIPROT"),
      row.names = all_prot_ids
    )

    return(anno_df)
  }

  if (is.na(annotation_db)) {
    return(anno_df)
  }
}


# add_annotation_hippie() --------

#' add_annotation_hippie ()
#'
#' @param anno_df annotation dataframe (for corresponding species in hippie)
#' @param ppi_hippie variable defined by ppis_hippie in get_networkdata_hippie()
#' @param species  from which species does the data come from
#'
#'
#'@return ppi_hippie with annotation columns for each interactor (for corresponding species in hippie)
#'
#'@export
#'
#'
#' @examples
#' \dontrun{
#'
#' db_hippie_df <- get_networkdata_hippie(species = "Homo_sapiens",
#'                                        version = "current",
#'                                        cache = TRUE,
#'                                        get_annotation = FALSE,
#'                                        add_annotation = FALSE
#'                                        )
#'
#' db_hippie_anno_df <- get_annotation_hippie(species = "Homo sapiens",
#'                                            version = "current"
#'                                            )
#'
#' db_hippie_ppi_anno_df <- add_annotation_hippie(ppi_hippie = db_hippie_df,
#'                                                anno_df = db_hippie_anno_df,
#'                                                species = "Homo_sapiens"
#'                                                )
#'
#'}

add_annotation_hippie <- function(ppi_hippie,
                                  anno_df,
                                  species
                                  ) {
  #adding Ensembl
  ppi_hippie$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_hippie$GeneSymbol_A, anno_df$genesymbol)]
  ppi_hippie$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_hippie$GeneSymbol_B, anno_df$genesymbol)]

  #adding Uniprot
  ppi_hippie$Uniprot_A <-
    anno_df$uniprot_id[match(ppi_hippie$GeneSymbol_A, anno_df$genesymbol)]
  ppi_hippie$Uniport_B <-
    anno_df$uniprot_id[match(ppi_hippie$GeneSymbol_B, anno_df$genesymbol)]

  return(ppi_hippie)

}


# outside of function ----------
# build_graph_hippie() -----

#' build_graph_hippie()
#'
#' @param graph_data ppi data from hippie
#' @param output_format selection of different graph functions that can be used
#' @param min_score_threshold select ppis that are "confident" depending on the scoretype/value
#'
#' @importFrom igraph graph.data.frame simplify
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#'
#' db_hippie_df <- get_networkdata_hippie(species = "Homo_Sapiens",
#'                                        version = "current"
#'                                       )
#'
#' db_hippie_graph <- build_graph_hippie(graph_data = db_hippie_df,
#'                                       output_format = "igraph",
#'                                       min_score_threshold = "0.35"
#'                                       )
#' db_hippie_graph
#' }
#'
#'

build_graph_hippie <- function (graph_data,
                                output_format = "igraph",
                                min_score_threshold = NULL ){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$`Confidence Value` <- as.numeric(graph_data$`Confidence Value`)

  # Erstelle das Histogramm mit 50 bins (breaks)
  hist(graph_data$`Confidence Value`, breaks = 50)

  #select ppi data >= minimal score
  if (!is.null(min_score_threshold)){
    graph_data_processed <- graph_data[graph_data$`Confidence Value` >= min_score_threshold, ]
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




