% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{reorder_list}
\alias{reorder_list}
\title{Helper function to reorder lists according to cvFolds}
\usage{
reorder_list(a_list, a_0, validRows, n_SL = 1, grn_ind = FALSE, n,
  for_se_cv = FALSE)
}
\arguments{
\item{a_list}{Structured list of nuisance parameters}

\item{a_0}{Treatment levels}

\item{validRows}{List of rows of data in validation data for
each split.}

\item{n_SL}{Number of super learners. If >1, then predictions
are averaged}

\item{grn_ind}{Structure of grn call is slightly different}

\item{n}{Sample size}

\item{for_se_cv}{Is this being used to average over
cross-validated standard errors? Affects index of \code{a_list}.}
}
\description{
Helper function to reorder lists according to cvFolds
}
