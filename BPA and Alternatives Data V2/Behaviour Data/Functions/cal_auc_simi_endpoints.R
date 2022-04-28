
#' Create expressions for matching
#'
#' @param segment a vector of items that need to be matched in the time_end
#' @param segment_name the name of the segment
#'
#' @return an expression
#' @keywords internal
#' @noRd
segment_2_expr <- function(segment, segment_name) {
  
  segment <- rlang::enquo(segment)
  segment_name <- rlang::enquo(segment_name)
  
  result <- rlang::expr(
    as.character(.data$time_end) %in% as.character(!!segment) ~ as.character(!!segment_name)
  )
  return(result)
}

cal_trapezoid <- function(x, y) {
  result <- c(NA, (y[-1] + y[-length(x)])*diff(x)/2)
  #result <- c(NA, (x$value[-1] + x$value[-nrow(x)])*diff(x$time_end)/2)
  return(result)
}

#' Add auc column on the tibble
#'
#' @param d 
#'
#' @return a tibble with an new auc column
#' @keywords internal
#' @noRd
add_trapezoid_value <- function(d) {
  result <- d %>%
    tidyr::nest(data = c("time_end", "value")) %>%
    dplyr::mutate(
      auc = purrr::map(data, ~ cal_trapezoid(.x$time_end, .x$value))
    ) %>%
    tidyr::unnest(cols = c("data", "auc"))
  return(result)
}

#' Add the endpoint column on the tibble
#'
#' Rows that do not have a match based on the expressions will be removed.
#' (endpoint is NA)
#'
#' @param d 
#' @param exprs 
#'
#' @return a tibble with an new endpoint column
#' @keywords internal
#' @noRd
#' 
add_endpoint_name <- function(d, exprs) {
  
  result <- d %>%
    dplyr::mutate(
      endpoint = dplyr::case_when(
        !!!exprs
      )
    ) %>%
    dplyr::filter(!is.na(.data$endpoint))
  return(result)
}

#' Summarize AUC values by endpoints
#'
#' @param d 
#'
#' @return
#' @keywords internal
#' @noRd
#'
#' 

summarize_auc_by_endpoint <- function(d) {
  result <- d %>%
    dplyr::group_by_at(dplyr::vars(-tidyselect::one_of(c("value", "time_end", "auc")))) %>%
    dplyr::summarize(endpoint_value = sum(.data$auc, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  return(result)
}

#' Check the dataset after excluding the time_end and value columns
#'
#' @param d 
#'
#' @return the original d if the data is okay
#' @keywords internal
#' @noRd
#' 
check_dataset <- function(d) {
  
  result <- d %>% 
    dplyr::distinct_at(dplyr::vars(-tidyselect::one_of(c("time_end", "value")))) %>% 
    dplyr::group_by_at(dplyr::vars(everything())) %>% 
    dplyr::count() %>%
    dplyr::filter(.data$n != 1) %>%
    dplyr::ungroup()
  
  if (nrow(result) > 0) {
    message(result)
    rlang::abort("id is not unique")
  }
  return(d)
}

#' Create endpoints based on area under the curve (AUC) in time segments
#' 
#' The endpoint value is the AUC in a time segment.
#' 
#' @param d a tibble with columns. Other columns can be included
#' \itemize{
#'   \item plate_id:
#'   \item embryo_id:
#'   \item time_end: numeric values
#'   \item value: the y axis (value > 0)
#' }
#' @param segments a named list. The name is the endpoint name. The value is the time_end that needs to be matched.
#'
#' @return a tibble with 
#' @export
#' 
#' 
create_auc_endpoints <- function(d, segments) {
  
  
  exps <- purrr::imap(segments, segment_2_expr)
  d <- check_dataset(d)
  
  d1 <- add_trapezoid_value(d)
  d2 <- add_endpoint_name(d1, exps)
  result <- summarize_auc_by_endpoint(d2)
  
  return(list(summary = result))
}

cosine <- function(x) {
  y <- t(x) %*% x
  res <- 1 - y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
  return(res)
}


#' Calculate similarity/distance matrix
#'
#' @param d a tibble with column names are the embryo_id
#' @param metric one of the four kinds
#'
#' @return a matrix that column/row names are the embryo_id
#' @keywords internal
#' @noRd
#' 
#' 
cal_simi_mat <- function(d, metric = c("pearson", "spearman", "euclidean", "cosine")) {
  mat <- as.matrix(d)
  
  if (metric == "pearson") {
    result <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  } else if (metric == "spearman") {
    result <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  } else if (metric == "euclidean") {
    result <- dist(t(mat), diag = TRUE, upper = TRUE, method = "euclidean")
  } else if (metric == "cosine" ) {
    result <- 1 - cosine(mat)
  }
  result <- as.matrix(result)
  return(result)
}

#' Add a simi_matrix column per plate_id + endpoint
#'
#' @param d 
#' @param metric 
#'
#' @return a tibble with simi_matrix + plate_id + endpoint
#' @keywords internal
#' @noRd
#' 
add_simi_matrix_col <- function(d, metric) {
  
  # plate_id + endpoint
  result <- d %>%
    dplyr::select(-.data$is_VC) %>%
    # make sure it is in the right sequence
    dplyr::arrange(.data$plate_id, .data$embryo_id, .data$endpoint, .data$time_end) %>% 
    tidyr::pivot_wider(., names_from = "embryo_id", values_from = "value") %>%
    dplyr::select(-.data$time_end) %>%
    tidyr::nest(data = -tidyselect::one_of("plate_id", "endpoint")) %>%
    dplyr::mutate(
      simi_matrix = purrr::map(data,  cal_simi_mat, metric = metric)
    ) %>%
    dplyr::select(-data)
  
  return(result)
}

#' Add a vc_ids column
#'
#' @param d 
#'
#' @return a tibble with plate_id + vc_ids, the vc_ids is a vector 
#' @keywords internal
#' @noRd
#' 
add_vc_ids_col <- function(d) {
  # plate_id
  result <- d %>%
    # is_VC = TRUE
    dplyr::filter(.data$is_VC) %>%
    dplyr::distinct(.data$plate_id, .data$embryo_id) %>%
    tidyr::nest(vc_ids = .data$embryo_id) %>%
    dplyr::mutate(vc_ids = purrr::map(.data$vc_ids, unlist, use.names = FALSE))
  return(result)
}

#' Collect the pairs that are related to the reference embryo_id
#' 
#' 1. exclude itself-itself
#' 2. include itself-other vc
#'
#' @param ... 
#'
#' @return a vector
#' @keywords internal
#' @noRd

collect_related_pair <- function(...) {
  l <- list(...)
  r_name <- l$embryo_id
  c_name <- l$vc_ids
  mat <- l$simi_matrix
  
  c_name <- c_name[!c_name %in% r_name]
  return(mat[r_name, c_name]) 
}

#' Add endpoint_value column that has the similarity distribution
#' 
#' endpoint + plate_id + embryo_id
#'
#' @param d 
#' @param matd 
#' @param VCd 
#'
#' @return
#' @keywords internal
#' @noRd

add_simi_distr_col <- function(d, matd, VCd) {
  #plate_id + endpoint + embryo_id
  result <- d %>%
    dplyr::select(-one_of(c("value", "time_end"))) %>% 
    dplyr::distinct() %>%
    dplyr::left_join(matd, by = c("plate_id", "endpoint")) %>%
    dplyr::left_join(VCd, by = "plate_id")
  
  result <- result %>%
    dplyr::mutate(endpoint_value = purrr::pmap(., collect_related_pair))
  return(result)
}


#' Calculate the mean of a vector with similarity values with some special treatments
#'
#' @param vec 
#'
#' @return a number
#' @keywords internal
#' @noRd
#' 
cal_mean_simi_vals <- function(vec) {
  out <- unlist(vec)
  out[is.na(out)] <- 0 # not moving at all; so that pearson/cosine/spearman = NA
  out[out < 0] <- 0 # for spearman & pearson
  result <- mean(out)
  return(result)
}

#' Summarize the similarity value per endpoint + plate_id + embryo_id
#'
#' @param d 
#'
#' @return
#' @keywords internal
#' @noRd
#' 
summarize_simi_by_endpoint <- function(d) {
  result <- d %>%
    dplyr::select(-.data$simi_matrix, -.data$vc_ids) %>%
    dplyr::mutate(
      endpoint_value  = purrr::map_dbl(.data$endpoint_value, cal_mean_simi_vals)
    )
  return(result)
}

unnest_simi_distri <- function(d) {
  result <- d %>% 
    dplyr::mutate(
      temp = purrr::map(.data$endpoint_value, enframe, name = "vs_embryo_id", value = "endpoint_value")) %>%
    dplyr::select(-.data$simi_matrix, -.data$vc_ids, -.data$endpoint_value) %>%
    tidyr::unnest(cols = .data$temp)
  return(result)
}

#' Create endpoints based on similarity of movement patterns in time segments
#'
#' @param d 
#' \itemize{
#'   \item plate_id:
#'   \item embryo_id:
#'   \item is_VC: vehicle control, it will be converted to logical
#'   \item time_end: numeric values
#'   \item value: the y axis (value > 0)
#' }
#' @param segments 
#' @param metric only one is allowed from pearson, spearman, euclidean, cosine
#'
#' @return a named list (summary + distribution) 
#' @export
#' 
#' 
create_simi_endpoints <- function(d, segments, metric = c("pearson", "spearman", "euclidean", "cosine")) {
  
  
  
  d <- check_dataset(d)
  d <- d %>% dplyr::mutate(is_VC = as.logical(.data$is_VC))
  
  # add endpoint
  exps <- purrr::imap(segments, segment_2_expr)
  d1 <- add_endpoint_name(d, exps)
  
  # plate_id + endpoint
  d_mat <- add_simi_matrix_col(d1, metric = metric)
  
  # plate_id
  d_ctrl <- add_vc_ids_col(d1)
  
  # plate_id + embryo_id + endpoint
  result1 <- add_simi_distr_col(d1, matd = d_mat, VCd = d_ctrl)
  
  # summarize
  result2 <- summarize_simi_by_endpoint(result1)
  result3 <- unnest_simi_distri(result1)
  
  return(list(summary = result2, distribution = result3))

}