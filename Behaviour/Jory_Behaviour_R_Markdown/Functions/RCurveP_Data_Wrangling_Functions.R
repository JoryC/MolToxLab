library(tidyverse)
library(DescTools)

median2 <- function(x){
  x <- na.omit(x)
  x1 <- sort(x)
  if (IsOdd(length(x1))) {
    return(median(x1))
  }
  if (all(x <= 0) & !IsOdd(length(x1))) {
    return(suppressWarnings(min(tail(x1, -floor(length(x1)/2)))))
  }
  if (!IsOdd(length(x1))) {
    return(suppressWarnings(max(head(x1, -floor(length(x1)/2)))))
  }
}


summarise_rcurvep_results <-
  function(resp_set,
           act_set,
           act_summary,
           reject_hit_conf_under = 0.5) {
    bmd_summary <- act_summary %>%
      group_by(chemical) %>%
      summarise(
        POD_med = POD_med,
        POD_cil = POD_cil,
        hit_confidence = hit_confidence,
        lowest_conc = lowest_conc
      ) # To merge with other sets
    
    result <- resp_set %>%
      inner_join(act_set) %>%
      group_by(chemical, conc, sample_id) %>%
      summarise(
        endpoint = endpoint,
        resp = resp,
        corrected_resp = corrected_resp,
        correction = correction,
        POD = POD,
        highest_conc = highest_conc,
        hit = hit,
        lowest_conc = lowest_conc
      ) %>%
      inner_join(bmd_summary) %>%
      ungroup() %>%
      group_by(chemical) %>%
      mutate(resp = if_else(correction == 1, true = corrected_resp, false = resp),
             lowest_conc = highest_conc-(number_of_dose_groups_per_chem-1)) %>%
      #dplyr::select(-corrected_resp) %>%
      mutate(median_POD = if_else(
        condition = POD == lowest_conc,
        true = unique(highest_conc),
        false = POD
      )) %>%
      mutate(median_POD = if_else(
        condition = is.na(median2(median_POD)),
        true = unique(highest_conc),
        false = median2(median_POD)
      )) %>%
      mutate(median_POD = if_else(
        condition = median_POD == Inf,
        true = unique(highest_conc),
        false = median2(median_POD)
      )) %>%
      mutate(
        median_POD = if_else(
          condition = hit_confidence < reject_hit_conf_under,
          true = unique(highest_conc),
          false = median_POD
        )
      ) %>%
      mutate(mean_POD = if_else(
        condition = POD == lowest_conc,
        true = unique(highest_conc),
        false = POD
      )) %>%
      mutate(mean_POD = if_else(
        condition = is.na(mean(mean_POD, na.rm = TRUE)),
        true = unique(highest_conc),
        false = mean(mean_POD, na.rm = TRUE)
      )) %>%
      mutate(
        mean_POD = if_else(
          condition = hit_confidence < reject_hit_conf_under,
          true = unique(highest_conc),
          false = mean_POD
        )
      ) %>%
      mutate(cil_POD = if_else(
        condition = POD == lowest_conc,
        true = unique(highest_conc),
        false = POD
      )) %>%
      mutate(cil_POD = if_else(
        condition = is.na(quantile(
          cil_POD, prob = c(0.05), na.rm = TRUE
        )),
        true = unique(highest_conc),
        false = quantile(cil_POD, prob = c(0.05), na.rm = TRUE)
      )) %>%
      mutate(
        cil_POD = if_else(
          condition = hit_confidence < reject_hit_conf_under,
          true = unique(highest_conc),
          false = cil_POD
        )
      ) %>%
      mutate(ciu_POD = if_else(
        condition = POD == lowest_conc,
        true = unique(highest_conc),
        false = POD
      )) %>%
      mutate(ciu_POD = if_else(
        condition = is.na(quantile(
          ciu_POD, prob = c(0.95), na.rm = TRUE
        )),
        true = unique(highest_conc),
        false = quantile(ciu_POD, prob = c(0.95), na.rm = TRUE)
      )) %>%
      mutate(
        ciu_POD = if_else(
          condition = hit_confidence < reject_hit_conf_under,
          true = unique(highest_conc),
          false = ciu_POD
        )
      )
    
    return(result)
  }



summarise_fitted_results <- function (resp_set,
                                      act_set,
                                      fit_set,
                                      act_summary,
                                      reject_hit_conf_under = 0.5) {
  temp <- resp_set %>%
    inner_join(act_set) %>%
    inner_join(fit_set) %>%
    group_by(chemical, conc, sample_id) %>%
    summarise(
      resp = resp,
      fitted_resp = fitted_resp,
      POD = POD,
      lowest_conc = lowest_conc,
      highest_conc = highest_conc,
      endpoint = endpoint,
      hit = hit,
      hill_fit = hill_fit,
      hill_aic = hill_aic,
      hill_tp = hill_tp,
      hill_ga = hill_ga,
      hill_gw = hill_gw
    ) %>%
    inner_join(act_summary) %>%
    dplyr::select(
      endpoint,
      chemical,
      sample_id,
      conc,
      resp,
      fitted_resp,
      hill_fit,
      hill_aic,
      hill_tp,
      hill_ga,
      hill_gw,
      lowest_conc,
      highest_conc,
      POD,
      POD_cil,
      POD_ciu,
      POD_med,
      hit,
      hit_confidence
    )
  
  result <- temp %>%
    group_by(chemical, conc, sample_id) %>%
    summarise(
      resp = unique(resp),
      fitted_resp = unique(fitted_resp),
      hill_fit = unique(hill_fit),
      hill_aic = unique(hill_aic),
      hill_tp = unique(hill_tp),
      hill_ga = unique(hill_ga),
      hill_gw = unique(hill_gw),
      lowest_conc = unique(lowest_conc),
      highest_conc = unique(highest_conc),
      POD = unique(POD),
      POD_cil = unique(POD_cil),
      POD_ciu = unique(POD_ciu),
      POD_med = unique(POD_med),
      hit = unique(hit),
      hit_confidence = unique(hit_confidence)
    ) %>%
    mutate(lowest_conc = highest_conc-(number_of_dose_groups_per_chem-1)) %>%
    mutate(median_POD = if_else(
      condition = POD == lowest_conc,
      true = unique(highest_conc),
      false = POD
    )) %>%
    mutate(median_POD = if_else(
      condition = median_POD >= highest_conc,
      true = unique(highest_conc),
      false = median_POD
    )) %>%
    mutate(median_POD = if_else(
      condition = is.na(median2(median_POD)),
      true = unique(highest_conc),
      false = median2(median_POD)
    )) %>%
    mutate(median_POD = if_else(
      condition = hit_confidence < reject_hit_conf_under,
      true = unique(highest_conc),
      false = median_POD
    )) %>%
    mutate(median_POD = if_else(
      condition = median_POD %in% c(NA, Inf),
      true = unique(highest_conc),
      false = median_POD
    )) %>%
    mutate(mean_POD = if_else(
      condition = POD == lowest_conc,
      true = unique(highest_conc),
      false = POD
    )) %>%
    mutate(mean_POD = if_else(
      condition = POD >= highest_conc,
      true = unique(highest_conc),
      false = mean_POD
    )) %>%
    mutate(mean_POD = if_else(
      condition = is.na(mean(mean_POD, na.rm = TRUE)),
      true = unique(highest_conc),
      false = mean(mean_POD, na.rm = TRUE)
    )) %>%
    mutate(mean_POD = if_else(
      condition = hit_confidence < reject_hit_conf_under,
      true = unique(highest_conc),
      false = mean_POD
    )) %>%
    mutate(cil_POD = if_else(
      condition = POD == lowest_conc,
      true = unique(highest_conc),
      false = POD
    )) %>%
    mutate(cil_POD = if_else(
      condition = POD >= highest_conc,
      true = unique(highest_conc),
      false = cil_POD
    )) %>%
    mutate(cil_POD = if_else(
      condition = is.na(quantile(
        cil_POD, prob = c(0.05), na.rm = TRUE
      )),
      true = unique(highest_conc),
      false = quantile(cil_POD, prob = c(0.05), na.rm = TRUE)
    )) %>%
    mutate(
      cil_POD = if_else(
        condition = hit_confidence < reject_hit_conf_under,
        true = unique(highest_conc),
        false = cil_POD
      )
    ) %>%
    mutate(ciu_POD = if_else(
      condition = POD == lowest_conc,
      true = unique(highest_conc),
      false = POD
    )) %>%
    mutate(ciu_POD = if_else(
      condition = POD >= highest_conc,
      true = unique(highest_conc),
      false = ciu_POD
    )) %>%
    mutate(ciu_POD = if_else(
      condition = is.na(quantile(
        ciu_POD, prob = c(0.95), na.rm = TRUE
      )),
      true = unique(highest_conc),
      false = quantile(ciu_POD, prob = c(0.95), na.rm = TRUE)
    )) %>%
    mutate(
      ciu_POD = if_else(
        condition = hit_confidence < reject_hit_conf_under,
        true = unique(highest_conc),
        false = ciu_POD
      )
    )
  
  return(result)
  
  
}



bmd_results <- function(summary_dat) {
  summary_dat %>%
    group_by(chemical) %>%
    summarise(
      median_POD = unique(median_POD),
      POD_med = unique(POD_med),
      mean_POD = unique(mean_POD),
      ciu_POD = unique(ciu_POD),
      cil_POD = unique(cil_POD),
      hit_confidence = unique(hit_confidence),
      lowest_conc = unique(lowest_conc),
      highest_conc = unique(highest_conc)
    ) %>%
    mutate(lowest_conc = highest_conc-(number_of_dose_groups_per_chem-1))
}


confident_hits <-
  function(summary_dat, reject_hit_conf_under = 0.5) {
    confident_hit_result <- summary_dat %>%
      filter(hit == 1, hit_confidence >= reject_hit_conf_under) %>%
      group_by(chemical, sample_id) %>%
      mutate(
        POD_diff_cil = abs(POD - cil_POD),
        POD_diff_ciu = abs(POD - ciu_POD),
        POD_diff_median = abs(POD - median_POD),
        POD_diff_med = abs(POD - POD_med),
        POD_diff_mean = abs(POD - mean_POD),
        smallest_diff = min(
          c(POD_diff_cil, POD_diff_ciu, POD_diff_mean, POD_diff_median, POD_diff_med)
        ),
        conf_POD = if_else(
          condition = smallest_diff == POD_diff_cil,
          true = "cil",
          false = if_else(
            condition = smallest_diff == POD_diff_ciu,
            true = "ciu",
            false = if_else(
              condition = smallest_diff == POD_diff_med,
              true = "med",
              false = if_else(
                condition = smallest_diff == POD_diff_median,
                true = "median",
                false = "mean"
              )  
            )
          )
        )
      ) %>%
      ungroup() %>%
      group_by(conf_POD) %>%
      group_split(.keep = TRUE)
    
    temp_name <- NULL
    for (i in 1:length(confident_hit_result)) {
      temp <- unique(confident_hit_result[[i]]$conf_POD)
      temp_name[i] <- temp
    }
    
    names(confident_hit_result) <- c(temp_name)
    
    return(confident_hit_result)
  }