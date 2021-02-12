#' Container for transcript annotation data
#'
#' Used to serve data to gene annotation tracks. Wraps \code{\link{GenomicRanges}} objects.
#' Annotation obtained from columns \code{Gene} (gene symbols) and \code{Exons} (exon start and end locations).
#'
#' @docType class
#' @seealso EpivizData
#' @seealso register,OrganismDb
#'
EpivizTxInfoData <- setRefClass("EpivizTxInfoData",
  contains="EpivizGeneInfoData"
)

.valid.EpivizTxInfoData.metadata <- function(x) {
  mdata <- mcols(x$.object)
  nms <- names(mdata)
  requiredNames <- c("Gene", "Transcript","Exons")
  if (any(!requiredNames %in% nms))
    return("'metadata' must contain columns 'Gene', 'Transcript' and 'Exons'")
  
  if (is(mdata$Gene, "Rle") && !is.character(runValue(mdata$Gene)))
    return("'Gene' must be a 'character' vector or Rle, or 'CharacterList'")
  if (!is(mdata$Gene, "Rle") && !(is.character(mdata$Gene) || is(mdata$Gene, "CharacterList")))
    return("'Gene' must be a 'character' vector or Rle, or 'CharacterList'")
  
  if (is(mdata$Transcript, "Rle") && !is.character(runValue(mdata$Transcript)))
    return("'Transcript' must be a 'character' vector or Rle, or 'CharacterList'")
  if (!is(mdata$Transcript, "Rle") && !(is.character(mdata$Transcript) || is(mdata$Transcript, "CharacterList")))
    return("'Transcript' must be a 'character' vector or Rle, or 'CharacterList'")
  
  if (!is(mdata$Exons, "IRangesList"))
    return("'Exons' must be an 'IRangesList'")
  
  NULL
}

EpivizTxInfoData$methods(
  get_default_chart_type = function() { "TranscriptTrack" },
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
                        transcript_id=as.character(.self$.object$Transcript[cur_hits]),
                        exon_starts=unname(lapply(start(.self$.object$Exons)[cur_hits], paste, collapse=",")),
                        exon_ends=unname(lapply(end(.self$.object$Exons)[cur_hits],paste,collapse=",")))
      out[[col]] <- cur_out
    }
    out
  },
  get_default_chart_type_html = function() {
    "epiviz-json-transcript-track"
  },
  get_metadata_columns = function() {
    c("transcript_id", "gene", "exon_starts","exon_ends")
  }
)