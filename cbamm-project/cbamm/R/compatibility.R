
#' Convert from metafor rma object
#' @export
from_metafor <- function(rma_obj) {
  data.frame(
    study_id = seq_len(length(rma_obj$yi)),
    yi = rma_obj$yi,
    se = sqrt(rma_obj$vi),
    vi = rma_obj$vi
  )
}

#' Convert from meta package object
#' @export  
from_meta <- function(meta_obj) {
  data.frame(
    study_id = meta_obj$studlab,
    yi = meta_obj$TE,
    se = meta_obj$seTE,
    vi = meta_obj$seTE^2
  )
}

#' Convert to metafor format
#' @export
to_metafor <- function(cbamm_data) {
  std_data <- standardize_cbamm_data(cbamm_data)
  list(yi = std_data$yi, vi = std_data$vi, sei = std_data$se)
}

#' Convert to meta format  
#' @export
to_meta <- function(cbamm_data) {
  std_data <- standardize_cbamm_data(cbamm_data)
  data.frame(
    TE = std_data$yi,
    seTE = std_data$se,
    studlab = std_data$study_id
  )
}

