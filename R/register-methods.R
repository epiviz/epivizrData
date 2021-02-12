#' Generic method to register data to the data server
#' 
#' @param object The object to register to data server
#' @param columns Name of columns containing data to register
#' @param ... Additonal arguments passed to object constructors
#' @return Object inheriting from \code{\link{EpivizData}} class
#' @export
#' 
#' @examples
#' library(GenomicRanges)
#' # create an example GRanges object
#' gr <- GRanges("chr10", IRanges(start=1:1000, width=100), score=rnorm(1000))
#' # this returns an EpivizData object without adding to data manager
#' # this is not the preferred way of creating these object, but is shown
#' # here for completeness.
#' ms_obj <- epivizrData:::register(gr, type="bp", columns="score") 
#' 
#' server <- epivizrServer::createServer(port=7123L)
#' data_mgr <- epivizrData::createMgr(server)
#'
#' # This adds a data object to the data manager
#' data_mgr$add_measurements(gr, "example_gr", type="bp", columns="score")
#'
setGeneric("register", signature=c("object"), 
           function(object, columns=NULL, ...) standardGeneric("register"))

setGeneric("reorderIfNeeded", signature=c("object"),
           function(object, ...) standardGeneric("reorderIfNeeded"))

setGeneric("coerceIfNeeded", signature=c("object"), 
           function(object, ...) standardGeneric("coerceIfNeeded"))

# TODO: add a sort check
setMethod("reorderIfNeeded", "GenomicRanges",
          function(object, ...) {
            stranded <- any(strand(object) != "*")
            if (stranded) {
              oobj <- object
              strand(object) <- "*"
            }
            if (!S4Vectors::isSorted(object)) {
              order <- order(object)
              if (stranded) {
                object <- oobj[order,]
              } else {
                object <- object[order,]
              }
            }
            return(object)
          })

setMethod("coerceIfNeeded", "GenomicRanges", 
          function(object, ...) {
            if (!is(object, "GNCList")) {
              newobject <- as(object, "GNCList")
              mcols(newobject) <- mcols(object)
              object <- newobject
            }
            object
          })

setMethod("reorderIfNeeded", "RangedSummarizedExperiment",
          function(object, ...) {
            gr <- rowRanges(object)
            stranded <- any(strand(gr) != "*")
            if (stranded) {
              ogr <- gr
              strand(gr) <- "*"
            }
            
            if (!S4Vectors::isSorted(gr)) {
              order <- order(gr)
              object <- object[order,]
            }
            object
          })

setMethod("coerceIfNeeded", "RangedSummarizedExperiment",
          function(object, ...) {
            if (!is(rowRanges(object), "GNCList")) {
              rowRanges(object) <- as(rowRanges(object), "GNCList")
            }
            object
          })

#' @describeIn register Register a \code{\link{GenomicRanges}} object
#' @import GenomicRanges
#' @param type Which type of data object to register for a \code{\link{GenomicRanges}} object. \code{block}: only region data, \code{bp} base-pair resolution quantitative data (see \code{columns} argument), \code{geneInfo} information about gene location.
setMethod("register", "GenomicRanges",
          function(object, columns, type=c("block","bp","tx_info","gene_info"), ...) {
            type <- match.arg(type)
            object <- reorderIfNeeded(object)
            object <- coerceIfNeeded(object)
            
            dev <- switch(type,
                          block = EpivizBlockData$new(object=object, ...),
                          bp = EpivizBpData$new(object=object, columns=columns, ...),
                          tx_info = EpivizTxInfoData$new(object=object, ...),
                          gene_info = EpivizGeneInfoData$new(object=object, ...)
            )
            return(dev)
          })

#' @describeIn register Register a \code{\link{RangedSummarizedExperiment}} object
#' @import SummarizedExperiment
#' @param assay Which assay in object to register
#' @param metadata Additional metadata about features
setMethod("register", "RangedSummarizedExperiment",
          function(object, columns = NULL, assay = 1, metadata = NULL) {
            object <- reorderIfNeeded(object)
            rowRanges(object) <- coerceIfNeeded(rowRanges(object))
            
            mcol_names <- names(mcols(rowRanges(object)))
            if (is.null(metadata) && !is.null(mcol_names)) {
              metadata <- mcol_names
            }
            if (!is.null(metadata) && any(!metadata %in% mcol_names)) {
              stop("invalid metadata")
            }
            EpivizFeatureData$new(object = object, 
                                  columns = columns, 
                                  assay = assay, 
                                  metadata = metadata)
          })

#' @describeIn register Register an \code{\link{ExpressionSet}} object
#' @param annotation Character string indicating platform annotation (only hgu133plus2 supported for now)
setMethod("register", "ExpressionSet",
          function(object, columns, annotation = NULL, assay="exprs") {
            if (is.null(annotation) || missing(annotation))
              annotation <- annotation(object)
            
            if (annotation != 'hgu133plus2') {
              stop("only 'hgu133plus2' affy chips supported for now")
            }
            
            # make GRanges object with appropriate info
            probeids <- featureNames(object)
            annoName <- paste0(annotation, ".db")
            
            if (!require(annoName, character.only=TRUE)) {
              stop("package '", annoName, "' is required")
            }
            
            res <- suppressWarnings(AnnotationDbi::select(get(annoName), keys=probeids, columns=c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND"), keytype="PROBEID"))
            dups <- duplicated(res$PROBEID)
            res <- res[!dups,]
            
            drop <- is.na(res$CHR) | is.na(res$CHRLOC) | is.na(res$CHRLOCEND)
            res <- res[!drop,]
            
            gr <- GRanges(seqnames = paste0("chr",res$CHR),
                          strand = ifelse(res$CHRLOC>0, "+","-"),
                          ranges = IRanges::IRanges(start=abs(res$CHRLOC), end=abs(res$CHRLOCEND)))
            
            mcols(gr)[,"SYMBOL"] <- res$SYMBOL
            mcols(gr)[,"PROBEID"] <- res$PROBEID
            
            mat <- assayDataElement(object, assay)[!drop,]
            if (missing(columns) || is.null(columns))
              columns <- colnames(mat)
            
            if (any(!(columns %in% colnames(mat))))
              stop("'columns' not found is 'assayDataElement(object, assay)'")
            
            mat <- mat[,columns]
            colnames(mat) <- columns
            
            if (!all(columns %in% rownames(pData(object)))) {
              pd <- data.frame(dummy=character(length(columns)))
              rownames(pd) <- columns
            } else {
              pd <- pData(object)[columns,]
            }
            sumexp <- SummarizedExperiment(assays = SimpleList(mat),
                                           rowRanges = gr,
                                           colData = DataFrame(pd))
            
            register(sumexp, 
                     columns = columns, 
                     assay = 1,
                     metadata = c("PROBEID","SYMBOL"))
          })

# 
# setMethod("register", "BigWigFile",
#           function(object, ...) {
#             dev <- EpivizWigData$new(file=object, ...)
#             return(dev)
# })
#           
# setMethod("register", "GAlignments",
#           function(object, coverage.only=TRUE, ...) {
#             if (!coverage.only) {
#               stop("'coverage.only' must be 'TRUE'. Only coverage supported for GAlignments for now.")
#             }
#             cov <- coverage(object)
#             register(as(cov,"GRanges"), columns="score", type="bp", ...)
# })
# 
# setMethod("register", "BamFile",
#           function(object, coverage.only=TRUE, ...) {
#             if (!coverage.only) {
#               stop("'coverage.only' muse be 'TRUE'. Only coverage supported for BamFiles for now.")
#             }
#             cov <- coverage(object)
#             register(as(cov,"GRanges"), columns="score", type="bp", ...)
# })

.register_txdb <- function(object,
                           kind = c("gene", "tx"),
                           keepSeqlevels=NULL, ...)
{
  message("creating gene annotation (it may take a bit)\n")
  kind <- match.arg(kind)
  gr <- .make_gene_info_gr(object, kind, keepSeqlevels)
  args <- list(...)
  if (!is.null(args$type)) {
    register(gr, ...)
  } else {
    register(gr, type = "gene_info", ...)
  }
}

#' @describeIn register Register an \code{\link{OrganismDb}} object
#' @import OrganismDbi
#' @param kind Make gene or transcript annotation (only gene supported for now)
#' @param keepSeqlevels character vector indicating seqlevels in object to keep   
setMethod("register", "OrganismDb", .register_txdb)

#' @describeIn register Register a \code{\link{TxDb}} object
#' @import GenomicFeatures
setMethod("register", "TxDb", .register_txdb)

#' @describeIn register Register an \code{\link{EnsDb}} object
#' @import ensembldb
setMethod("register", "EnsDb", .register_txdb)

# setMethod("register", "BEDFile", function(object, ...) {
#   gr <- import.bed(object)
#   register(gr, type="block", ...)
# })

#' @describeIn register Register an \code{\link{data.frame}}
setMethod("register", "data.frame",
          function(object, columns, ...) {
            if(!("chr" %in% colnames(object))) {
              object$chr <- rep("dataframe", times=nrow(object))
            }
            if(!("start" %in% colnames(object))) {
              object$start <- 1:nrow(object)
            }
            if(!("end" %in% colnames(object))) {
              object$end <- object$start
            }
            rownames(object) <- c()
            
            gdf <- makeGRangesFromDataFrame(object, keep.extra.columns = TRUE)
            register(gdf, 
                     type = "bp",
                     columns = columns, ...)
          })
