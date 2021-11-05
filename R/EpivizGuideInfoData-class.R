#' Container for guide annotation data
#'
#' Used to serve data to transcript annotation tracks. Wraps \code{\link{GenomicRanges}} objects.
#'
#' @docType class
#' @seealso EpivizData
#'
EpivizGuideInfoData <- setRefClass("EpivizGuideInfoData",
  contains="EpivizGeneInfoData"
)

.valid.EpivizGuideInfoData.metadata <- function(x) {
  mdata <- mcols(x$.object)
  nms <- names(mdata)
  requiredNames <- c("ID", "pam_site", "spacer_20mer", "pam", "cut_site", "nuclease")
  if (any(!requiredNames %in% nms))
    return("'metadata' must contain columns 'ID', 'pam_site', 'spacer_20mer', 'pam', 'cut_site', 'nuclease'")
  
  if (is(mdata$ID, "Rle") && !is.character(runValue(mdata$ID)))
    return("'ID' must be a 'character' vector or Rle, or 'CharacterList'")
  
  if (!is(mdata$Exons, "IRangesList"))
    return("'Exons' must be an 'IRangesList'")
  
  NULL
}


EpivizGuideInfoData$methods(
  get_default_chart_type = function() { "GuideTrack" },
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
                        ID=as.character(.self$.object$ID[cur_hits]),
                        pam_site=as.character(.self$.object$pam_site[cur_hits]),
                        spacer_20mer=as.character(.self$.object$spacer_20mer[cur_hits]),
                        pam=as.character(.self$.object$pam[cur_hits]),
                        cut_site=as.character(.self$.object$cut_site[cur_hits]),
                        nuclease=as.character(.self$.object$nuclease[cur_hits])
                        )
      out[[col]] <- cur_out
    }
    out
  },
  get_default_chart_type_html = function() {
    "epiviz-json-guide-track"
  },
  get_metadata_columns = function() {
    c("ID", "pam_site", "spacer_20mer","pam", "cut_site", "nuclease")
  }
)