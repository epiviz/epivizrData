#' Class providing data manager for epiviz app
#'
#' @import methods IRanges
#' @importClassesFrom epivizrServer EpivizServer
#' @importClassesFrom GenomeInfoDb Seqinfo
#' @export 
EpivizDataMgr <- setRefClass("EpivizDataMgr",
  fields = list(
    .ms_list = "environment",
    .ms_idCounter = "integer",
    .server = "EpivizServer",
    .seqinfo = "Seqinfo"
  ),
  methods=list(
    initialize=function(server=epivizrServer::createServer(), ...) {
      .self$.server <- server
      .self$.ms_list <- new.env(parent=emptyenv())
      .self$.ms_idCounter <- 0L
    },
    num_datasources = function() { length(ls(.self$.ms_list)) },
    show = function() {
      "Print manager information to screen"
      cat("Epiviz data manager object:\n")
      cat("Server: ")
      .self$.server$show(); cat("\n")
 
#       st <- .self$list_measurements()
#       if (length(st)>0) {
#         cat("Measurements:\n")
#         print(st); cat("\n")
    },
    is_server_closed = function() {
      "Check if underlying server is closed, <logical>"
      is.null(.self$.server) || .self$.server$is_closed()
    }
  )
)

# measurement management methods
EpivizDataMgr$methods(
  .clear_datasourceGroup_cache = function(ms_obj) {
    if(!is(ms_obj, "EpivizData")) {
      stop("'ms_obj' must be an 'EpivizData' object")
    }
    if (!exists(ms_obj$get_id(), envir = .self$.ms_list, inherits = FALSE))
      stop("did not find object")
    
    if (!.self$is_server_closed()) {
      callback <- function(response_data) {
        cat(ms_obj$get_id(), " datasourceGroup caches cleared\n")
        invisible()
      }
      
      request_data <- list(action = "clearDatasourceGroupCache",
                           datasourceGroup = ms_obj$get_id())
      .self$.server$send_request(request_data, callback)
    }
  },
  is_ms_connected = function(ms_obj_or_id) {
    "check if measurement object was properly added to JS app"
    ms_object <- .self$.get_ms_object(ms_obj_or_id)
    ms_id <- ms_object$get_id()
    isTRUE(.self$.ms_list[[ms_id]]$connected)
  },
  add_measurements = function(obj, datasource_name=NULL, datasource_origin_name=deparse(substitute(obj)), send_request = TRUE, ...) {
    "register measurements in data manager"
    if (missing(datasource_name) || is.null(datasource_name)) {
      datasource_name <- datasource_origin_name
    }
    
    if (!is.character(datasource_name)) {
      stop("data source name has to be a string: ", datasource_name)
    }
    ms_object <- epivizrData:::register(obj, ...)

    .self$.ms_idCounter <- .self$.ms_idCounter + 1L
    ms_id <- sprintf("%s_%d", datasource_name, .self$.ms_idCounter)
    ms_object$set_id(ms_id)
    ms_object$set_name(datasource_name)
    ms_object$set_source_name(datasource_obj_name)
    ms_object$set_mgr(.self)
    
    measurements <- ms_object$get_measurements()
    
    ms_record <- list(measurements=measurements,
                      name=datasource_name, 
                      source_name=datasource_origin_name,
                      obj=ms_object, 
                      connected=FALSE)
    assign(ms_id, ms_record, envir=.self$.ms_list)
    
    send_request <- !.self$is_server_closed() && isTRUE(send_request)
    if (send_request) {
      .callback <- function(response_data) {
        .self$.ms_list[[ms_id]]$connected <- TRUE
        cat("Measurement ", datasource_name, " added to application and connected\n")
      }
      request_data <- list(action="addMeasurements",
                           measurements=epivizrServer::json_writer(lapply(measurements, as.list)))
      .self$.server$send_request(request_data, .callback)
    }
    ms_object
  },
  get_measurements = function() {
    out <- .emptyEpivizMeasurement()

    measurements <- list()
    ids <- ls(.self$.ms_list)
    if (length(ids)>0) {
      for (id in ids) {
        ms_record <- .self$.ms_list[[id]]
        ms <- ms_record$obj$get_measurements()
        for (cur_ms in ms) {
          out <- .appendEpivizMeasurement(out, cur_ms)
        }
      }
    }

    as.list(out)
  },
  .get_ms_object=function(ms_obj_or_id) {
    ms_obj <- NULL
    if (is.character(ms_obj_or_id)) {
      # passed the id instead of the object
      id <- ms_obj_or_id
      if (!exists(id, envir=.self$.ms_list, inherits=FALSE)) {
        stop("measurement with id ", id, " not found")
      }
      ms_obj <- .self$.ms_list[[id]]$obj
    } else {
      ms_obj <- ms_obj_or_id
    }
    ms_obj
  },
  rm_measurements=function(ms_obj_or_id) {
    "remove registered measurments from a given data object"
    ms_obj <- .get_ms_object(ms_obj_or_id)
    
    if (!is(ms_obj, "EpivizData")) {
      stop("'ms_obj' must be an 'EpivizData' object")
    }
    
    id <- ms_obj$get_id()
    if (!exists(id, envir=.self$.ms_list, inherits=FALSE)) {
      stop("measurement with id ", id, " not found")
    }
    
    ms_record <- .self$.ms_list[[id]]
    ms_name <- ms_record$name
    ms <- ms_record$obj$get_measurements()
    rm(list=id, envir=.self$.ms_list)
    
    if(isTRUE(ms_record$connected)) {
      callback <- function(response_data) {
        cat("measurement object ", ms_name, " removed and disconnected\n")
      }
      request_data <- list(action="removeMeasurements",
                           measurements=epivizrServer::json_writer(lapply(ms, as.list)))
      .self$.server$send_request(request_data, callback)
    }
    invisible()
  },
  rm_all_measurements = function() {
    "remove all registered measurements"
    ids <- ls(.self$.ms_list)
    if (length(ids)>0) {
      for (id in ids) {
        .self$rm_measurements(id)
      }
    }
  },
  update_measurements = function(ms_obj_or_id, new_object, send_request = TRUE) {
    "update the underlying data object for a registered measurement (given by object or id)"
    ms_obj <- .self$.get_ms_object(ms_obj_or_id)
    if (!is(ms_obj, "EpivizData"))
      stop("ms_obj must be of class 'EpivizData'")
    ms_obj$update(new_object, send_request=send_request)
    invisible()
    
  },
  list_measurements = function() {
    "make a printable list of registered measurements"
    if (.self$num_datasources() == 0) {
      return(data.frame())
    }
    
    ids <- ls(.self$.ms_list)
    nms <- sapply(ids, function(id) .self$.ms_list[[id]]$name)
    lens <- sapply(ids, function(id) length(.self$.ms_list[[id]]$obj$.object))
    connected <- ifelse(sapply(ids, function(id) .self$.ms_list[[id]]$connected), "*", "")
    columns <- sapply(ids, function(id) paste0(.self$.ms_list[[id]]$obj$.columns, collapse=","))
    
    data.frame(id=ids,
               name=nms,
               length=lens,
               connected=connected,
               columns=columns,
               stringsAsFactors=FALSE, row.names=NULL)
  }
)

# data fetch methods
EpivizDataMgr$methods(
  .find_datasource = function(datasource) {
    if (!exists(datasource, .self$.ms_list, inherits=FALSE)) {
      stop("cannot find datasource", datasource)
    }
    ms_obj <- .self$.ms_list[[datasource]]$obj
  },
  get_rows = function(chr, start, end, metadata, datasource) {
    if (is.null(chr) || is.null(start) || is.null(end)) {
      query <- NULL
    } else {
      query <- GRanges(chr, ranges=IRanges(start, end))
    }
    ms_obj <- .self$.find_datasource(datasource)
    ms_obj$get_rows(query, metadata)
  },
  get_values = function(chr, start, end, datasource, measurement) {
    if (is.null(chr) || is.null(start) || is.null(end)) {
      query <- NULL
    } else {
      query <- GRanges(chr, ranges=IRanges(start, end))
    }
    ms_obj <- .self$.find_datasource(datasource)
    ms_obj$get_values(query, measurement)
  }
)

# seqinfo methods
EpivizDataMgr$methods(
  add_seqinfo = function(seqinfo, keep_seqlevels=NULL, send_request=TRUE) {
    if(!is(seqinfo, "Seqinfo")) {
      stop("'seqinfo' must be of class 'Seqinfo'")
    }
    if (!is.null(keep_seqlevels)) {
      seqinfo <- keepSeqlevels(seqinfo, keep_seqlevels)
    }
    
    .self$.seqinfo <- seqinfo

    send_request <- !.self$is_server_closed() && isTRUE(send_request)
    if (send_request) {
      callback <- function(response_data) {
        cat("Seqinfos added to application\n")
      }
      seqlengths <- seqlengths(.self$.seqinfo)+1
      seqinfo_list <- mapply(function(seqname, seqlength) list(seqname, 1, seqlength), 
                             names(seqlengths), seqlengths,
                             SIMPLIFY=FALSE, USE.NAMES=FALSE)
      request_data <- list(action="addSeqInfos",
                           seqInfos=epivizrServer::json_writer(seqinfo_list))
      .self$.server$send_request(request_data, callback)
    }
    invisible()
  },
  rm_seqinfo = function(send_request=TRUE) {
    if (length(.self$.seqinfo) == 0) {
      return(invisible())
    }
    
    seqnames <- seqnames(.self$.seqinfo)
    .self$.seqinfo <- Seqinfo()
    
    send_request <- !.self$is_server_closed() && isTRUE(send_request)
    if (send_request) {
      callback <- function(response_data) {
        cat("Seqinfos deleted from app\n")
      }
      request_data <- list(action="removeSeqNames",
                           seqNames=epivizrServer::json_writer(seqnames))
      .self$.server$send_request(request_data, callback)
    }
    invisible()
  },
  get_seqinfo = function() {
    seqlengths <- seqlengths(.self$.seqinfo)+1
    seqinfo_list <- structure(lapply(seqlengths, function(x) c(1,x)), names=names(seqlengths))
    seqinfo_list
  }
)