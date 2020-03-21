#' Methylation reference for cell types in blood
#'
#' Provides expected methylation levels for cell type
#' markers in blood for estimating cell type fractions
#'
#' @docType data
#'
#' @usage data(blood.ref.333cpgs)
#'
#' @format An object of class `"matrix"` with 333 CpG sites
#'         as rows and cell types as columns
#'
#' @keywords datasets
#'
#' @references Teschendorff et al. 2017, BMC Bioinformatics.
#'
#' @examples
#' data(blood.ref.333cpgs)
#' cell.types <- colnames(blood.ref.333cpgs)
"blood.ref.333cpgs"


#' Koestler Methylation Dataset
#'
#' Includes over 400k methylation sites from blood for 18
#' individuals with FACS-based cell fractions for major
#' blood cell types available. 
#'
#' @docType data
#'
#' @usage data(koestler)
#'
#' @format An object of class `"list"`. Slot \code{X} contains
#'         methylation values in \code{matrix} format. Slots
#'         \code{C1} and \code{C2} are \code{NULL} as no
#'         covariates were used in these analyses.
#'         Slot \code{W.facs} contains FACS-based cell fractions
#'         in a \code{data.frame} format.
#'
#' @keywords datasets
#'
#' @references Koestler DC, Jones MJ, Usset J, Christensen BC, Butler RA, Kobor MS, Wiencke JK, Kelsey KT: Improving cell mixture deconvolution by identifying optimal DNA methylation libraries (IDOL). BMC Bioinformatics 2016, 17:120.
#'
#' @examples
#' data(koestler)
#' samples <- colnames(koestler$X)
"koestler"
