% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/register-methods.R
\name{register}
\alias{register}
\alias{register,GenomicRanges-method}
\alias{register,RangedSummarizedExperiment-method}
\alias{register,ExpressionSet-method}
\alias{register,OrganismDb-method}
\alias{register,TxDb-method}
\alias{register,EnsDb-method}
\alias{register,data.frame-method}
\title{Generic method to register data to the data server}
\usage{
register(object, columns = NULL, ...)

\S4method{register}{GenomicRanges}(object, columns, type = c("block", "bp", "tx_info", "gene_info"), ...)

\S4method{register}{RangedSummarizedExperiment}(object, columns = NULL, assay = 1, metadata = NULL)

\S4method{register}{ExpressionSet}(object, columns, annotation = NULL, assay = "exprs")

\S4method{register}{OrganismDb}(object, kind = c("gene", "tx"), keepSeqlevels = NULL, ...)

\S4method{register}{TxDb}(object, kind = c("gene", "tx"), keepSeqlevels = NULL, ...)

\S4method{register}{EnsDb}(object, kind = c("gene", "tx"), keepSeqlevels = NULL, ...)

\S4method{register}{data.frame}(object, columns = NULL, ...)
}
\arguments{
\item{object}{The object to register to data server}

\item{columns}{Name of columns containing data to register}

\item{...}{Additonal arguments passed to object constructors}

\item{type}{Which type of data object to register for a \code{\link{GenomicRanges}} object. \code{block}: only region data, \code{bp} base-pair resolution quantitative data (see \code{columns} argument), \code{geneInfo} information about gene location.}

\item{assay}{Which assay in object to register}

\item{metadata}{Additional metadata about features}

\item{annotation}{Character string indicating platform annotation (only hgu133plus2 supported for now)}

\item{kind}{Make gene or transcript annotation (only gene supported for now)}

\item{keepSeqlevels}{character vector indicating seqlevels in object to keep}
}
\value{
Object inheriting from \code{\link{EpivizData}} class
}
\description{
Generic method to register data to the data server
}
\section{Methods (by class)}{
\itemize{
\item \code{GenomicRanges}: Register a \code{\link{GenomicRanges}} object

\item \code{RangedSummarizedExperiment}: Register a \code{\link{RangedSummarizedExperiment}} object

\item \code{ExpressionSet}: Register an \code{\link{ExpressionSet}} object

\item \code{OrganismDb}: Register an \code{\link{OrganismDb}} object

\item \code{TxDb}: Register a \code{\link{TxDb}} object

\item \code{EnsDb}: Register an \code{\link{EnsDb}} object

\item \code{data.frame}: Register an \code{\link{data.frame}}
}}

\examples{
library(GenomicRanges)
# create an example GRanges object
gr <- GRanges("chr10", IRanges(start=1:1000, width=100), score=rnorm(1000))
# this returns an EpivizData object without adding to data manager
# this is not the preferred way of creating these object, but is shown
# here for completeness.
ms_obj <- epivizrData:::register(gr, type="bp", columns="score") 

server <- epivizrServer::createServer(port=7123L)
data_mgr <- epivizrData::createMgr(server)

# This adds a data object to the data manager
data_mgr$add_measurements(gr, "example_gr", type="bp", columns="score")

}
