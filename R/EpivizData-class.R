setClassUnion("EpivizDataMgrOrNULL", c("EpivizDataMgr", "NULL"))

#' Data container for epiviz data server
#' 
#' @import methods
#' @import S4Vectors
#' @export
EpivizData <- setRefClass("EpivizData",
  contains="VIRTUAL",
  fields=list(
    .object="ANY",
    .mgr="EpivizDataMgrOrNULL",
    .id="character",
    .name="character",
    .source_name="character",
    .columns="ANY",
    .ylim="ANY",
    .cur_query="ANY",
    .cur_hits="ANY"
  ),
  methods=list(
    initialize = function(object = GNCList(GRanges()), 
                          columns=NULL, 
                          ylim=NULL, ...) {
      .self$.object <- object
      .self$.columns <- columns

      if (!.self$.check_columns(columns))
        stop("Invalid 'columns' argument")

      if (is.null(.self$.columns))
        .self$.columns <- .self$.infer_columns()

      naIndex <- .self$.infer_nas()
      if (length(naIndex)>0) {
        .self$.object <- .self$.object[-naIndex,]
      }
      
      if (!is.null(ylim)) {
        if (!.self$.check_limits(ylim))
          stop("invalid 'ylim' argument")
        .self$.ylim <- ylim
      } else {
        .self$.ylim <- .self$.infer_limits()
      }

      .self$.cur_query <- NULL
      .self$.cur_hits <- NULL
      .self$.mgr <- NULL
      callSuper(...)
    },
    .infer_nas = function() {
      integer()
    },
    .check_columns = function(columns) {
      is.null(columns)
    },
    .infer_columns = function() {
      NULL
    },
    .check_limits = function(ylim) {
      is.null(ylim)
    },
    .infer_limits = function() {
      NULL
    },
    .check_class = function(new_object) {
      stop(".check_class called on VIRTUAL object")
    },
    update = function(new_object, send_request=TRUE) {
      "Update underlying data object with new object"
      if (!.self$.check_class(new_object)) {
        stop("class of 'new_object' is not equal to class of current 'object'")
      }

      original_object <- .self$.object
      .self$.object <- new_object
      
      if (!is.null(.self$.columns)) {
        if (!.self$.check_columns(.self$.columns)) {
          .self$.object <- original_object
          stop("columns not found in 'new_object'")
        }
        .self$.ylim <- .self$.infer_limits()
      }

      .self$.object <- reorderIfNeeded(.self$.object)
      .self$.object <- coerceIfNeeded(.self$.object)
      
      na_index <- .self$.infer_nas()
      if (length(na_index) > 0) {
        .self$.object <- .self$.object[-na_index,]
      }
      
      if (send_request && !is.null(.self$.mgr))
        .self$.mgr$.clear_datasourceGroup_cache(.self)
      invisible()
    },
    # TODO: use accessor functions for this
    get_id = function() {
      "Get id provided by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.id
    },
    set_id = function(id) {
      "Set id, used by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.id <- id
      invisible()
    },
    
    # TODO: use accessor functions for this
    get_name = function() { 
      "Get datasource name, usually set by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.name 
    },
    set_name=function(name) {
      "Set datasource name, usually set by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.name <- name
      invisible()
    },
    
    # TODO: use accessor functions for this
    get_source_name = function() {
      "Get original datasource name provided by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.source_name
    },
    set_source_name = function(source_name) {
      "Set original datasource name, used by manager \\code{\\link{EpivizDataMgr-class}}"
      .self$.source_name <- source_name
      invisible()
    },
    
    set_limits = function(ylim) {
      "Set plotting limits for continuous data"
      if (!.self$.check_limits(ylim))
          stop("'invalid' limits argument")
      .self$.ylim <- ylim
    }, 
    get_measurements = function() {
      "Get description of measurements served by this object"
      stop("'get_measurements' called on virtual class object")
    },
    get_default_chart_type = function() {
      "Get name of default chart type for this data type"
      stop("'get_default_chart_type' called on virtual class object")
    },
    parse_measurement = function(ms_id=NULL) {
      "Parse a measurement description for data served by this object"
      stop("'parse_measurement' called on virtual class object")
    },
    set_mgr = function(mgr) {
      "Set data manager, \\code{\\link{EpivizDataMgr-class}}"
      if (!is(mgr, "EpivizDataMgr"))
        stop("'mgr' must be of class 'EpivizDataMgr'")
      
      .self$.mgr <- mgr
      invisible()
    },
    show = function() {
      "Print information about this object"
      cat(class(.self), "object", .self$.id, "\n")
      methods::show(.self$.object)
      cat("\n\tcolumns:", paste(.self$.columns, collapse=","),"\n")
      cat("\tlimits:\n")
      print(.self$.ylim)
    }
  )
)

#####
# validity
.valid.EpivizData.columns <- function(x) {
  if(!x$.check_columns(x$.columns))
    return("invalid 'columns' slot")
  NULL
}

.valid.EpivizData <- function(x) {
  c(.valid.EpivizData.columns(x))
}

S4Vectors::setValidity2("EpivizData", .valid.EpivizData)

#######
# get data
EpivizData$methods(
  .get_hits = function(query) {
    if (!is(query, "GRanges"))
      stop("'query' must be a GRanges object")
    if (length(query) != 1) {
      stop("'query' must be of length 1")
    }

    if (is.null(.self$.cur_query) || 
        !identical(unname(query), unname(.self$.cur_query))) 
    {
      .self$.cur_query <- query
      olaps <- suppressWarnings(GenomicRanges::findOverlaps(query, .self$.object, select="all"))
      .self$.cur_hits <- subjectHits(olaps)

      if (length(.self$.cur_hits) == 0) {
        return(invisible())
      }

      if (!S4Vectors::isSorted(start(.self$.object)[.self$.cur_hits])) {
        stop("these should be ordered by now...")
      }
      # this is here to make insertion on JS cache work
      .self$.cur_hits <- seq(min(.self$.cur_hits), max(.self$.cur_hits))
    }
    invisible()
  },
  get_rows = function(query, metadata, useOffset = FALSE) {
    "Get genomic interval information overlapping query <\\code{\\link{GenomicRanges}}> region"
    if (is.null(query)) {
      out <- list(globalStartIndex=NULL, useOffset=FALSE,
                  values=list(id=list(),
                              start=list(),
                              end=list(),
                              metadata=.self$.get_metadata(integer(), metadata)))
      return(out)
    }

    .self$.get_hits(query)
    if (length(.self$.cur_hits) == 0) {
      out <- list(globalStartIndex=NULL, useOffset=FALSE,
                  values=list(id=list(),
                    start=list(),
                    end=list(),
                    metadata=.self$.get_metadata(.self$.cur_hits, metadata)))
    } else {
      if (!useOffset) {
        out <- list(globalStartIndex=.self$.cur_hits[1],
                  useOffset=FALSE,
                  values=list(
                    id=.self$.cur_hits,
                    start=start(.self$.object)[.self$.cur_hits],
                    end=end(.self$.object)[.self$.cur_hits],
                    metadata=.self$.get_metadata(.self$.cur_hits, metadata)
                   )
                )
      } else {
        st <- start(.self$.object)[.self$.cur_hits]
        stDiff <- diff(st)
        end <- end(.self$.object)[.self$.cur_hits]
        endDiff <- diff(end)

        out <- list(globalStartIndex=.self$.cur_hits[1],
                    useOffset=TRUE,
                    values=list(
                      id=.self$.cur_hits,
                      start=c(st[1], stDiff),
                      end=c(end[1], endDiff),
                      metadata=.self$.get_metadata(.self$.cur_hits, metadata)
                     )
                )
        }
    }
    # make sure we are wrapping lists correctly
    if (length(out$values) > 0 && length(out$values$id) == 1) {
      for (slotName in names(out$values)) {
        if (slotName != "metadata")
          out$values[[slotName]] <- list(out$values[[slotName]])
      }
    }
    return(out)
  },
  .get_values_from_hits = function(curHits, measurement, round) {
    numeric()
  },
  get_values=function(query, measurement, round=TRUE) {
    "Get measurement values for features overlapping query region <\\code{\\link{GenomicRanges}}"
    if (is.null(query)) {
      out <- list(globalstartIndex=NULL, values=list())
      return(out)
    }

    .self$.get_hits(query)
    if (length(.self$.cur_hits) == 0) {
      out <- list(globalStartIndex=NULL, values=list())
    } else {
      out <- list(globalStartIndex=.self$.cur_hits[1],
                  values=.self$.get_values_from_hits(.self$.cur_hits, measurement, round=round))
      # check we are wrapping lists correctly
      if (length(out$values) == 1) {
        out$values <- list(out$values)
      }
    }
    return(out)
  }
)



