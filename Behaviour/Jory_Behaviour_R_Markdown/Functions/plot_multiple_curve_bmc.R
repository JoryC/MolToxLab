library(here)
library(scales)
library(tidyverse)
library(gridExtra)

#' Plot multiple chemical concentration response data annotated with BMC. 
#' six chemical and its responses on x endpoints.
#' By using the gridExtra::marrageGrob & ggplot2::ggsave(), 
#' output curves can be saved in mult-page pdf format
#'
#' @param ldata a list of tibbles with columns: plot_id, endpoint, concs, resps, and batch
#' each tibble data has unique plot_id and details of the columns are explained as followed
#' \itemize{
#'   \item plot_id: the title of each small plot, usually the chemical ID
#'   \item endpoint: the individual panel on each small plot, the # of endpoints = # of panels
#'   \item concs: DMSO (concs = 0) can be included, in microM
#'   \item resps: (normalized) responses with baseline = 0
#'   \item batch: batch id, for example, one chemical could be tested on three different plates
#'   \item n_per_page: default = 6
#' }
#' @param bmc_data a tibble with columns: plot_id, endpoint, BMC, score, batch_active_fraction
#' It is assumed that there is 1-to-1 relationship between endpoint and active BMC. 
#' For the siuation of 1-2 (endpoint-BMC), which could happen for the chemical 
#' that have both increasing and decreasing effect (U-shape). 
#' It it suggested to summerize the BMC (e.g., most potent BMC as the representative) before running the script
#' \itemize{
#'   \item plot_id: same as above
#'   \item endpoint: same as above
#'   \item BMC (in microM); inactive has NA values
#'   \item score: hit_confidence score from Rcurvep [0-1] or NA
#'   \item flag: flag (e.g., assay interference or fitting issues) for annotation or NA
#'   \item batch_active_fraction: the fraction of active if tested in multiple batches [0-1]
#'   If < 0.5, and "&" is added next to the BMC
#' }
#' @param type a character string to indicate whether to use percent style (e.g., percent of mortality) or continuous style
#'
#' @return the output from gridExtra::marrangeGrob with 6 chemicals on a page
#' @export
#'
#' @examples
#' source(here("plot_multiple_curve_bmc.R"))
#' httr::GET("https://github.com/moggces/sample_files/blob/master/biobide_zebrafish_DT_resp_bmc_data.rds?raw=true", httr::write_disk(here("sample_data.rds"), overwrite =  TRUE))
#' #httr::GET("https://github.com/moggces/sample_files/blob/master/biobide_zebrafish_NT_resp_bmc_data.rds?raw=true", httr::write_disk(here("sample_data.rds"), overwrite =  TRUE))
#' sample_data <- readRDS(here("sample_data.rds"))
#' names(sample_data)
#' outputc <- curve_bmc_annotation_plot_wrapper(sample_data[["resp_data"]], sample_data[["bmc_data"]], type = "percent")
#' #for biobide_zebrafish_NT_resp_bmc_data.rds, use height = 10
#' ggsave(here("DT_curves.pdf"), outputc, width = 15, height = 7)
#' 
#' 
curve_bmc_annotation_plot_wrapper <- function(ldata, bmc_data, type = c("percent", "continuous"), n_per_page = 6) {
  
  # find the most common response direction in the dataset
  #direction <- get_resp_direction(ldata)
  
  # create the annotation on the plots  
  anno_text <- purrr::map(ldata, get_annotation_text, bmc_data = bmc_data)
  
  # generate plots
  mpl <- purrr::map2(ldata, anno_text, single_curve_bmc_annotation_plot, type = type)
  
  # arrange the plots 
  result <- gridExtra::marrangeGrob(mpl, nrow = 1, ncol = n_per_page, top = "")
  
  return(result)
}

get_resp_direction <- function(dd) {
  #mean_resp <- ldata %>% dplyr::bind_rows() %>% dplyr::pull(resps) %>% mean(., na.rm = TRUE)
  mean_resp <- dd %>% dplyr::pull(resps) %>% mean(., na.rm = TRUE)
  direction <- 1
  
  if (mean_resp < 0 & !is.na(mean_resp)) {
    direction <- '-1'
  }
  return(direction)
}

get_annotation_text <- function(dd, bmc_data) {
  
  direction <- get_resp_direction(dd)
  pseudo_resp <- Inf
  if (direction == -1) {
    pseudo_resp <- -Inf
  }
  
  result <- bmc_data %>%
    dplyr::semi_join(
      dd, by = c("plot_id", "endpoint")
    ) %>%
    dplyr::mutate(
      #BMC = ifelse(score < 0.5, NA, BMC),
      score = ifelse(score < 0.5, NA, score)) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        !is.na(BMC) ~ stringr::str_c("BMC: ", scientific(BMC), "\nscore: ", stringr::str_replace_na(score)),
        TRUE ~ as.character(NA)
      ),
      label = dplyr::case_when(
        !is.na(BMC) & !is.na(flag) ~ str_c(label, "\nflag: ", stringr::str_replace_na(flag)),
        TRUE ~ label
      ),
      label = dplyr::case_when(
        !is.na(label) & batch_active_fraction < 0.5 ~ stringr::str_c(label, " &"),
        TRUE ~ label
      ) 
    ) %>%
    mutate(concs = -Inf, resps = pseudo_resp, batch = as.factor(1))
  return(result)
}

single_curve_bmc_annotation_plot <- function(dd, annotation, type) {
  
  
  if (type == "percent") {
    
    p <- get_percent_plot(dd, annotation)  
    
  } else if (type == "continuous") {
    p <- get_continuous_plot(dd, annotation)
  }
  return(p)
}

str_wrap_30 <- function(x) {
  stringr::str_wrap(x, 30)
}

get_percent_plot <- function(dd, annotation) {
  
  direction <- get_resp_direction(dd)
  sel_concs <- dd %>% dplyr::pull(concs) %>% unique() %>% sort()
  sel_name <- dd %>% dplyr::pull(plot_id) %>% unique()
  
  # make sure batch as factor
  dd <- dd %>% dplyr::mutate(batch = as.factor(batch))
  
  # replace 0 as DMSO
  x_axis_label <- scales::scientific(sel_concs)
  if (0 %in% sel_concs) {
    x_axis_label <- c("DMSO", scales::scientific(sel_concs[-1]))
  }
  
  # increasing/decreasing resps
  limits <- c(0, 100)
  if (direction == -1) {
    limits <- c(-100, 0)
  }
  
  # position of the annotation
  position <- 1.5
  if (direction == -1) {
    position <- -0.1
  }
  
  #make concentration as factor
  dd <- dd %>% 
    mutate(concs = forcats::fct_reorder(as.factor(concs), concs)) 
  
  # dot + line plot
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = concs, y = resps, color = batch, group = batch))
  p <- p + 
    ggplot2::geom_point(size = 3, alpha = 0.7) + 
    ggplot2::geom_line(alpha = 0.7) + 
    ggplot2::geom_text( data = annotation, 
                        mapping = ggplot2::aes(x = concs, y = resps, label = label),
                        hjust   = -0.1,
                        vjust   = position, 
                        color = "blue", size = 3 ) + 
    ggplot2::scale_y_continuous("response", limits = limits) + 
    ggplot2::scale_x_discrete(expression(paste(mu, "M", sep = "")), 
                              labels = x_axis_label) + 
    ggplot2::ggtitle(stringr::str_wrap(sel_name, width = 20)) + 
    ggplot2::guides(color = FALSE) + 
    ggplot2::theme_bw() + 
    ggplot2::theme(
      plot.title = ggplot2::element_text( size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.8), angle = 45, hjust = 1, vjust = 1)
    ) + 
    ggplot2::facet_wrap(~ endpoint, ncol = 1,  
                        labeller = labeller(endpoint = str_wrap_30)) 
}

get_continuous_plot <- function(dd, annotation) {
  
  direction <- get_resp_direction(dd)
  sel_concs <- dd %>% dplyr::pull(concs) %>% unique() %>% sort()
  sel_name <- dd %>% dplyr::pull(plot_id) %>% unique
  
  dd <- dd %>% dplyr::mutate(batch = as.factor(batch))
  
  
  x_axis_label <- scales::scientific(sel_concs)
  if (0 %in% sel_concs) {
    x_axis_label <- c("DMSO", scales::scientific(sel_concs[-1]))
  }
  
  # position of the annotation
  position <- 1.5
  if (direction == -1) {
    position <- -0.1
  }
  
  # boxplot for the batch (1st layer)
  # dot (2nd layer)
  # boxplot for the whole (3rd layer)
  dd <- dd %>% 
    dplyr::mutate(concs = forcats::fct_reorder(as.factor(concs), concs)) 
  
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = concs, y = resps, color = batch))
  p <- p + 
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5,  alpha = 0.7) + 
    ggplot2::geom_jitter(size = 1, alpha = 0.7) + 
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.5,  color = "black", alpha = 0.5) + 
    ggplot2::geom_text( data = annotation, 
                        mapping = ggplot2::aes(x = concs, y = resps, label = label),
                        hjust   = -0.1,
                        vjust   = position, 
                        color = "blue", size = 3 ) +
    ggplot2::scale_shape_manual(values = c("1" = 1, "0"  = 16)) + 
    ggplot2::scale_y_continuous("response") + 
    ggplot2::theme_bw() + 
    ggplot2::guides(shape = FALSE, color = FALSE) + 
    ggplot2::scale_x_discrete(expression(paste(mu, "M", sep = "")), 
                              labels = x_axis_label) + 
    ggplot2::ggtitle(stringr::str_wrap(sel_name, width = 20)) + 
    ggplot2::theme(
      plot.title = element_text( size = ggplot2::rel(1)),
      axis.text.x = element_text(size = ggplot2::rel(0.8), angle = 45, hjust = 1, vjust = 1)
    ) + 
    ggplot2::facet_wrap(~ endpoint, ncol = 1, scales = "free", 
                        labeller = labeller(endpoint = str_wrap_30))
  return(p)
}