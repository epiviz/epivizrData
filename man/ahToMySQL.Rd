% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ahToMySQL}
\alias{ahToMySQL}
\title{Utility function to import data to a MySQL database from Annotation Hub}
\usage{
ahToMySQL(ah, annotations = list(), ...)
}
\arguments{
\item{ah}{[AnnotationHub()] object with records to add to
database.}

\item{annotations}{A named list of lists (key/value pairs). Keys must 
be the AH ID for the corresponding record and the value is a named list
representing an annotation. An annotation is automatically inferred by
the record's metadata. Any annotation that is passed for a particular
record is concatenated to its inferred annotation. If the annotation has
a subtype column, it is used to name the data object being added to the db, 
otherwise the record's tags is used.}

\item{...}{arguments for toMySQL (connection, db_name, batch, index)}
}
\description{
Utility function to import data to a MySQL database from Annotation Hub
}
\examples{
\dontrun{
library(epivizrData)
library(AnnotationHub)
library(DBI)
library(RMySQL)

ah <- AnnotationHub()
db_annotations <- list()

# Query Patterns
roadmap <- "EpigenomeRoadMap"
bisulphite <- "bisulphite"

esophagus <- query(ah, c("esophagus", "roadmap", "bisulphite"))
eso_anno <- list(tissue="Digestive", subtype="Esophagus")
eso_id <- names(esophagus)
db_annotations[[eso_id]] <- eso_anno

connection <- dbConnect(MySQL(), host=host, user=user, password=pass)
db_name="my_database"

ahToMySQL(ah=record, annotations=db_annotations,
  connection=connection, db_name=db_name)
}
}
