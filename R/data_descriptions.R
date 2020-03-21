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
