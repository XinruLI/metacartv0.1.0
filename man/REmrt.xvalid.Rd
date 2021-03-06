% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SupportLookAhead.R
\name{REmrt.xvalid}
\alias{REmrt.xvalid}
\title{A function to perform cross-validation for RE meta-tree}
\usage{
REmrt.xvalid(mf, maxL, n.fold, minbucket, minsplit, delQ, lookahead)
}
\arguments{
\item{mf}{the data.frame to grow the tree}

\item{maxL}{the maximum number of splits}

\item{n.fold}{the number of folds to perform the cross-validation}

\item{minbucket}{the minimum number of the studies in a terminal node}

\item{minsplit}{the minimal number of studies in a parent node to be split}

\item{delQ}{the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.}

\item{lookahead}{an argument indicating whether to apply the "look-ahead" strategy when fitting the tree}
}
\value{
a cp table
}
\description{
A function to perform cross-validation for RE meta-tree
}
\keyword{internal}
