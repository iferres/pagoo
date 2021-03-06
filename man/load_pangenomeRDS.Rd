% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_functions.R
\name{load_pangenomeRDS}
\alias{load_pangenomeRDS}
\title{Load A Pagoo Pangenome}
\usage{
load_pangenomeRDS(file, pkg, ...)
}
\arguments{
\item{file}{The path to the pangenome `.RDS` file.}

\item{pkg}{The package to use to load the object. Shouldn't be necessary to
provide, but may be useful in some cases.}

\item{...}{Arguments to be passed to the pagoo object. \code{sep} and
\code{core_level} overwrite the values stored in the file.}
}
\value{
A \code{PgR6MS} class object, or a \code{PgR6M} object (with or
without sequences, respectively).
}
\description{
This function loads a pagoo pangenome from a `.RDS` file generated
by the `save_pangenomeRDS` method. Objects loaded by this functions keep their
states, i.e : dropped/recovered organisms are conserved, as well as the
`core_level` setted when the object was originally saved.
}
