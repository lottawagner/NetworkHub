

# get_networkdata_reactome() -----------

#' get_networkdata_reactome()
#'
#' @param species  default value = "taxid:9606(Homo sapiens)" - from which species does the data come from
#' @param version version of the data files in reactome
#' @param cache default value set to TRUE (automatically checks if the data file is already stored in the cache)
#' @param get_annotation creation of an annotation dataframe using AnnotationDbi packages, default value set to TRUE
#' @param add_annotation adding annotation to ppi dataframe, default value set to TRUE
#' @param ... 	further arguments passed to or from other methods
#'
#' @return ppis_reactome
#'
#' @importFrom vroom vroom
#' @export
#'
#' @examples
#' \dontrun{
#' db_reactome_df <- get_networkdata_reactome( species = "taxid:9606(Homo sapiens)",
#'                                             version = "current",
#'                                             cache = TRUE,
#'                                             get_annotation = FALSE,
#'                                             add_annotation = FALSE
#'                                            )
#' }

get_networkdata_reactome <- function(species = "taxid:9606(Homo sapiens)",
                                     version = "current",
                                     cache = TRUE,
                                     get_annotation = TRUE,
                                     add_annotation = TRUE,
                                     ...) {

  # list species is actualized for version reactome "current" (sept 2024)
  # UPDATEVERSION

  # check that the value for species is listed in reactome

  if (!(species %in% list_species_reactome)) { # if species is not in the list
    stop("Species not found as specified by reactome,",
         "please check some valid entries by running `list_species_reactome`") # stop function and print
  }


  # buildup of the resource location for the version and all
  ## elegantly done in another smaller utility function

  rname <- paste0(
    "reactome_v_",
    version
    )
  #UPDATEVERSION

  if (cache) {
    # tries to fetch from the cache
    message("Trying to fetch from cache...")
    network_file <- fetch_NetworkHub(rname)
  } # here we check if the file is already in the cache and use the function fetch_NetworkHub to take the corresponding file

  if (!cache | is.null(network_file)) {
    # retrieves the file for the first time
    message("Downloading to cache...")
    # buildup from "base" reactome url
    reactome_url <-
      urlmaker_reactome(
        version = version
      ) # if there is no entry for the corresponding file in the cache, we create the url using urlmaker_reactome

    # and cache_NetworkHub to cache the file from the url source
    network_file <- cache_NetworkHub(
      rname = rname,
      fpath = reactome_url
    )
  }

  # read in the resource, whether cached or freshly downloaded
  ppis_reactome <- vroom::vroom(network_file)
  #ppis_reactome <- head(read.delim(network_file, sep = " "))

  # If you want to check all species contained in the file
  #species_A <- unique(ppis_reactome$`Taxid interactor A`)
  #species_B <- unique(ppis_reactome$`Taxid interactor B`)
  #info_species_reactome <- union(species_A, species_B)


  #Filter for species

  ppis_reactome_filtered <- ppis_reactome[(ppis_reactome$`Taxid interactor A` == species ) &
                                          (ppis_reactome$`Taxid interactor B` == species ),]

  # rename columns

  #Uniprot
  colnames(ppis_reactome_filtered)[colnames(ppis_reactome_filtered) == "#ID(s) interactor A"] <- "Uniprot_A"
  colnames(ppis_reactome_filtered)[colnames(ppis_reactome_filtered) == "ID(s) interactor B"] <- "Uniprot_B"

  #extract Uniprot_id
  ppis_reactome_filtered$Uniprot_A <- str_extract(ppis_reactome_filtered$Uniprot_A, "uniprotkb:([A-Z0-9]+)")
  ppis_reactome_filtered$Uniprot_A <- gsub("uniprotkb:", "", ppis_reactome_filtered$Uniprot_A)

  ppis_reactome_filtered$Uniprot_B <- str_extract(ppis_reactome_filtered$Uniprot_B, "uniprotkb:([A-Z0-9]+)")
  ppis_reactome_filtered$Uniprot_B <- gsub("uniprotkb:", "", ppis_reactome_filtered$Uniprot_B)

  #Confidence score
  ppis_reactome_filtered$`Confidence value(s)` <- str_extract(ppis_reactome_filtered$`Confidence value(s)`, "reactome-score:([0-9\\.]+)")
  ppis_reactome_filtered$`Confidence value(s)`<- gsub("reactome-score:", "", ppis_reactome_filtered$`Confidence value(s)`)



  if (get_annotation && !is.na(annotation_db)){

    db_reactome_anno_df <- get_annotation_reactome( ppi_reactome = ppis_reactome_filtered,
                                                    species = species,
                                                    version = version
                                                  )

    message("...created annotation dataframe")

    if (add_annotation) {

      db_reactome_ppi_anno_df <- add_annotation_reactome( anno_df = db_reactome_anno_df,
                                                          ppi_reactome = ppis_reactome_filtered,
                                                          species = species
                                                        )

      message("...added annotation from *db_reactome_anno_df* to *db_reactome_ppi_anno_df*")

      return(db_reactome_ppi_anno_df)
    }

    if (!add_annotation){
      return(db_reactome_anno_df)
    }
  }

  if (!get_annotation) {
    if (add_annotation){
      stop("get_annotation must be = TRUE in order to add_annotation")
    }
  }
  return(ppis_reactome_filtered)
}


# Outside function ----------


#SPECIES NetworkHub uses from reactome
list_species_reactome <- c("taxid:9913(Bos taurus)",
                           "taxid:6239(Caenorhabditis elegans)",
                           "taxid:9615(Canis familiaris)",
                           "taxid:7227(Drosophila melanogaster)",
                           "taxid:562(Escherichia coli)",
                           "taxid:9031(Gallus gallus)",
                           "taxid:9606(Homo sapiens)",
                           "taxid:10090(Mus musculus)",
                           "taxid:10116(Rattus norvegicus",
                           "taxid:4932(Saccharomyces cerevisiae)",
                           "taxid:8355(Xenopus laevis)"
                          )


list_db_annotationdbi_reactome <- c("org.Bt.eg.db",
                                   "org.Ce.eg.db",
                                   "org.Cf.eg.db",
                                   "org.Dm.eg.db",
                                   "org.EcK12.eg.db",
                                   "org.Gg.eg.db",
                                   "org.Hs.eg.db",
                                   "org.Mm.eg.db",
                                   "org.Rn.eg.db",
                                   "org.Sc.sgd.db",
                                   "org.Xl.eg.db"
                                   )


reactome_db_annotations <- data.frame(species = list_species_reactome,
                                 anno_db_reactome = list_db_annotationdbi_reactome,
                                 row.names = list_species_reactome
)


# get_annotation_reactome() --------

#' get_annotation_reactome ()
#'
#' @param species  from which species does the data come from
#' @param version version of the data files in reactome
#' @param ppi_reactome variable defined by ppis_reactome in get_networkdata_reactome()
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats na.omit
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
#' @import org.Xl.eg.db
#'
#'@return ppi_reactome
#'
#'@export
#'
#'
#' @examples
#' \dontrun{
#' db_reactome_df <- get_networkdata_reactome( species = "taxid:9606(Homo sapiens)",
#'                                             version = "current",
#'                                             cache = TRUE,
#'                                             get_annotation = FALSE,
#'                                             add_annotation = FALSE
#'                                            )
#'
#' db_reactome_anno_df <- get_annotation_reactome( ppi_reactome = db_reactome_df,
#'                                                 species = "taxid:9606(Homo sapiens)"
#'                                                 )
#' }

get_annotation_reactome <- function(ppi_reactome,
                                    species,
                                    version) {

  # find database on corresponding species

  if (!(species %in% list_species_reactome)) { # if species is not in the list
    stop("Species not found as specified by reactome,",
         "please check some valid entries by running `list_species_reactome`") # stop function and print
  }

  annotation_db <-
    reactome_db_annotations$anno_db_reactome[match(species, reactome_db_annotations$species)]

  if (!is.na(annotation_db)) {

  all_prot_ids <- unique(na.omit(c(ppi_reactome$Uniprot_A, ppi_reactome$Uniprot_B)))

  anno_df <- data.frame(uniprot_id = all_prot_ids,
                        gene_symbol <- mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "UNIPROT",
                                              column = "SYMBOL"
                                              ),
                        ensembl_id <- mapIds(get(annotation_db),
                                              keys = all_prot_ids,
                                              keytype = "UNIPROT",
                                              column = "ENSEMBL"
                                            ),
                        entrez_ids <- mapIds(get(annotation_db),
                                             keys = all_prot_ids,
                                             keytype = "UNIPROT",
                                             column = "ENTREZID"
                                            ),
                        row.names = all_prot_ids
                        )

  return(anno_df)
  }

  if (is.na(annotation_db)) {
    message("Annotation database for the species is not implemented yet.\n",
            "Next time define add_annotation in get_networkdata_reactome(..., add_annotation = FALSE, ...)\n",
            "You will get ppis_reactome containing annotation for Uniprot_ and GeneSymbol_.")
    return(ppi_reactome)
  }
}


# add_annotation_reactome() --------
#' add_annotation_reactome ()
#'
#' @param anno_df annotation dataframe (for corresponding species in reactome)
#' @param ppi_reactome variable defined by ppis_reactome in get_networkdata_reactome()
#' @param species  from which species does the data come from
#'
#'
#' @return ppi_reactome with annotation columns for each interactor (for corresponding species in reactome)
#'
#' @export
#'
#'
#' @examples
#' \dontrun{
#' db_reactome_df <- get_networkdata_reactome( species = "taxid:9606(Homo sapiens)",
#'                                             version = "current",
#'                                             cache = TRUE,
#'                                             get_annotation = FALSE,
#'                                             add_annotation = FALSE
#'                                            )
#'
#' db_reactome_anno_df <- get_annotation_reactome( ppi_reactome = db_reactome_df,
#'                                                 species = "taxid:9606(Homo sapiens)"
#'                                                 )
#'
#' db_reactome_ppi_anno_df <- add_annotation_reactome( ppi_reactome = db_reactome_df,
#'                                                     anno_df = db_reactome_anno_df,
#'                                                     species = "taxid:9606(Homo sapiens)"
#'                                                     )
#'
#'}

add_annotation_reactome <- function(ppi_reactome,
                                    anno_df,
                                    species) {

  #adding GeneSymbol
  ppi_reactome$GeneSymbol_A <-
    anno_df$gene_symbol[match(ppi_reactome$Uniprot_A, anno_df$uniprot_id)]
  ppi_reactome$GeneSymbol_B <-
    anno_df$gene_symbol[match(ppi_reactome$Uniprot_B, anno_df$uniprot_id)]

  #adding Ensembl
  ppi_reactome$Ensembl_A <-
    anno_df$ensembl_id[match(ppi_reactome$Uniprot_A, anno_df$uniprot_id)]
  ppi_reactome$Ensembl_B <-
    anno_df$ensembl_id[match(ppi_reactome$Uniprot_B, anno_df$uniprot_id)]

  #adding Entrez
  ppi_reactome$Entrez_A <-
    anno_df$entrez_id[match(ppi_reactome$Uniprot_A, anno_df$uniprot_id)]
  ppi_reactome$Entrez_B <-
    anno_df$entrez_id[match(ppi_reactome$Uniprot_B, anno_df$uniprot_id)]

  return(ppi_reactome)
}




# build_graph_reactome() -----

#' build_graph_reactome()
#'
#' @param graph_data ppi data from reactome
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
#' db_reactome_df <- get_networkdata_reactome(species = "taxid:9606(Homo sapiens)",
#'                                            version = "current"
#'                                            )
#'
#' db_reactome_graph <- build_graph_reactome(graph_data = db_reactome_df,
#'                                         output_format = "igraph",
#'                                         min_score_threshold = "0.5")
#' db_reactome_graph #list of 5416
#' }
#'
#'


build_graph_reactome <- function (graph_data,
                              output_format = "igraph",
                              min_score_threshold = NULL ){

  #check on the clumns in your ppi data file
  colnames(graph_data)

  graph_data$`Confidence value(s)` <- as.numeric(graph_data$`Confidence value(s)`)

  # Erstelle das Histogramm mit 50 bins (breaks)
  hist(graph_data$`Confidence value(s)`, breaks = 50)

  #select ppi data >= minimal score
  if (!is.null(min_score_threshold)){
    graph_data_processed <- graph_data[graph_data$`Confidence value(s)` >= min_score_threshold, ]
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








