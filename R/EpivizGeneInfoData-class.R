#' Container for gene annotation data
#'
#' Used to serve data to gene annotation tracks. Wraps \code{\link{GenomicRanges}} objects.
#' Annotation obtained from columns \code{Gene} (gene symbols) and \code{Exons} (exon start and end locations).
#'
#' @docType class
#' @seealso EpivizData
#' @seealso register,OrganismDb
#'
EpivizGeneInfoData <- setRefClass("EpivizGeneInfoData",
  contains="EpivizTrackData",
  methods=list(
    initialize=function(...) {
      callSuper(...)
      .self$.columns <- NULL
    }
  )
)

.valid.EpivizGeneInfoData.ylim <- function(x) {
  if (!is.null(x$.ylim))
    return("'ylim' must be 'NULL'")
  NULL
}

.valid.EpivizGeneInfoData.metadata <- function(x) {
  mdata <- mcols(x$.object)
  nms <- names(mdata)
  requiredNames <- c("Gene","Exons")
  if (any(!requiredNames %in% nms))
    return("'metadata' must contain columns 'Gene' and 'Exons'")

  if (is(mdata$Gene, "Rle") && !is.character(runValue(mdata$Gene)))
    return("'Gene' must be a 'character' vector or Rle, or 'CharacterList'")
  if (!is(mdata$Gene, "Rle") && !(is.character(mdata$Gene) || is(mdata$Gene, "CharacterList")))
    return("'Gene' must be a 'character' vector or Rle, or 'CharacterList'")

  if (!is(mdata$Exons, "IRangesList"))
    return("'Exons' must be an 'IRangesList'")

  NULL
}

.valid.EpivizGeneInfoData <- function(x) {
  c(.valid.EpivizGeneInfoData.ylim(x),
    .valid.EpivizGeneInfoData.metadata(x))
}

S4Vectors::setValidity2("EpivizGeneInfoData", .valid.EpivizGeneInfoData)

EpivizGeneInfoData$methods(
  get_default_chart_type = function() { "GenesTrack" },
  get_measurements=function() {
    out <- list(EpivizMeasurement(
      id = .self$.id,
      name = .self$.name,
      type = "range",
      datasourceId = .self$.id,
      datasourceGroup = .self$.id,
      datasourceName = .self$.source_name,
      defaultChartType = .self$get_default_chart_type(),
      metadata=.self$get_metadata_columns()))
    out
  },
  get_rows=function(query, metadata) {
    out <- callSuper(query, metadata)
    if (length(.self$.cur_hits) == 0) {
      return(out)
    }

    out$values$strand <- as.character(strand(.self$.object)[.self$.cur_hits])
    out
  },
  .get_metadata=function(cur_hits, cur_metadata) {
    if (length(cur_hits) == 0) {
      out <- lapply(cur_metadata, function(x) list())
      names(out) <- cur_metadata
      return(out)
    }
    out <- vector("list", length(cur_metadata))
    names(out) <- cur_metadata
    for (col in cur_metadata) {
      cur_out <- switch(col,
                       gene=as.character(.self$.object$Gene[cur_hits]),
                       exon_starts=unname(lapply(start(.self$.object$Exons)[cur_hits], paste, collapse=",")),
                       exon_ends=unname(lapply(end(.self$.object$Exons)[cur_hits],paste,collapse=",")))
      out[[col]] <- cur_out
    }
    out
  },
  get_default_chart_type_html = function() {
    "epiviz-json-genes-track"
  },
  .get_col_data = function(query) {
    NULL
  },
  .get_sql_index_table_info = function(annotation) {
    if (is.null(annotation)) {
      annotation <- "NULL"
    }
    list(index_table="gene_data_index",
      values=list(paste0(
        "'", .self$get_id(), "'", ",", # measurement_id
        "'", .self$get_name(), "'", ",", # measurement_name
        "'", .self$get_id(), "'", ",", # location
        "'", .self$get_id(), "'", ",", # column_name
        0, ",", # min
        0, ",", # max
        "'", annotation, "'")
        )
      )
  },
  get_metadata_columns = function() {
    c("gene", "exon_starts","exon_ends")
  }
)
