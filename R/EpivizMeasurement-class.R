#' Class encapsulating a measurement description for epiviz app.
#' 
#' @exportClass EpivizMeasurement
EpivizMeasurement <- setClass("EpivizMeasurement",
  contains = "SparseEpivizMeasurement",
  slots = c(
    name = "character",
    type = "character",
    datasourceGroup = "character",
    defaultChartType = "character",
    dataprovider = "character",
    annotation = "ANY",
    minValue = "numeric",
    maxValue = "numeric",
    metadata = "ANY"
  ),
  prototype = prototype(
    name = character(),
    type = character(),
    datasourceGroup = character(),
    defaultChartType = character(),
    dataprovider = character(),
    annotation = NULL,
    minValue = as.numeric(NA),
    maxValue = as.numeric(NA),
    metadata = NULL
  )
)

#' Convert \code{\link{EpivizMeasurement}} object to \code{list}
#' 
#' @param x \code{\link{EpivizMeasurement}} object to coerce.
#' @return a \code{list} describing measurement object
#' @export
setMethod("as.list", signature(x="EpivizMeasurement"),
  function(x) {
    nms <- slotNames("EpivizMeasurement")
    out <- lapply(nms, function(slot_name) slot(x, slot_name))
    names(out) <- nms
    out
  }          
)

#' Display measurement datasourceId and id
#' 
#' @param object a \code{\link{EpivizMeasurement}} to display
#' @return A string describing measurement
#' @export
setMethod("show", signature(object="EpivizMeasurement"),
  function(object) {
    cat(paste0(object@datasourceId, ":", object@id))    
  }
)

#' Create empty Epiviz Measurement
#' @export
.emptyEpivizMeasurement <- function() {
  EpivizMeasurement(id=character(),
                    name=character(),
                    type=character(),
                    datasourceId=character(),
                    datasourceGroup=character(),
                    datasourceName=character(),
                    defaultChartType=character(),
                    dataprovider=character(),
                    annotation=list(),
                    minValue=numeric(),
                    maxValue=numeric(),
                    metadata=list())
}

setMethod(".appendEpivizMeasurement", "EpivizMeasurement",
          function(a, b) {
            nms <- slotNames("EpivizMeasurement")
            for (nm in nms) {
              cur_val <- slot(b, nm)
              if (is.list(slot(a, nm))) {
                cur_val <- list(cur_val)
              }
              if (!is.null(slot(b, nm))) {
                slot(a, nm) <- c(slot(a, nm), cur_val)
              } else {
                slot(a, nm) <- c(slot(a, nm), list(NULL))
              }
            }
            a
          })

setMethod(".serializeEpivizMeasurement", "EpivizMeasurement",
          function(x) {
            if (length(x@id) == 1) {
              nms <- slotNames("EpivizMeasurement")
              for (nm in nms) {
                slot(x, nm) <- list(slot(x, nm))
              }
            }
            
            epivizrServer::json_writer(as.list(x))
          })