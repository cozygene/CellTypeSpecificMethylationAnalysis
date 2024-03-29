% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_descriptions.R
\docType{data}
\name{koestler}
\alias{koestler}
\title{Koestler Methylation Dataset}
\format{An object of class `"list"`. Slot \code{X} contains
        methylation values in \code{matrix} format. Slots
        \code{C1} and \code{C2} are \code{NULL} as no
        covariates were used in these analyses.
        Slot \code{W.facs} contains FACS-based cell fractions
        in a \code{data.frame} format.}
\usage{
data(koestler)
}
\description{
Includes over 400k methylation sites from blood for 18
individuals with FACS-based cell fractions for major
blood cell types available.
}
\examples{
data(koestler)
samples <- colnames(koestler$X)
}
\references{
Koestler DC, Jones MJ, Usset J, Christensen BC, Butler RA, Kobor MS, Wiencke JK, Kelsey KT: Improving cell mixture deconvolution by identifying optimal DNA methylation libraries (IDOL). BMC Bioinformatics 2016, 17:120.
}
\keyword{datasets}
