### General line plot
# X is one of the parameters
# Y is reproducibility, sensitivity or specificity
# Lines represent the different measures of reproducibility (you can specify a subset)
# Facetted in grid by two of the parameters

### Extra settings
# For facetting or coloring or X axis, parameters can be:
#   bias, power, interlab variation, repro.power, etc
# Additionally, you can pass a filter to the functions, to be applied before
# All these setting are registered in the caption of the figure and in the filename

filtered_data = function (D, x, to_plot, facetting, types, to_include, aggregate_prev = F, which_agg_prev) {
  # browser()
  # Filter the data if a filter is given
  if (!rlang::quo_is_null(to_include)) {
    D = D %>% filter(rlang::eval_tidy(to_include, D))
  }
  
  # Filter the types to plot
  if (to_plot %in% c("reproducibility", "root mean square error"))
  { to_plot_suffix = "ReproRate" }
  else if (aggregate_prev) {
    if (to_plot == "sensitivity") { to_plot_suffix = c("TP", "FN") }
    else if (to_plot == "specificity") { to_plot_suffix = c("TN", "FP") }
  } else {
    if (to_plot == "sensitivity") { to_plot_suffix = "SENS" }
    else if (to_plot == "specificity") { to_plot_suffix = "SPEC" }
  }
  
  to_plot_cols = paste0(sort(rep(types, length(to_plot_suffix))), "_", to_plot_suffix)
  
  keep = c(to_plot_cols, facetting, x)
  if (aggregate_prev) keep = c(keep, which_agg_prev)
  
  D = D %>% select(all_of(keep))
  
  if (aggregate_prev & to_plot != "root mean square error") {
    agg_cols = c(facetting, x, which_agg_prev)
    D = D %>% group_by(across(all_of(agg_cols))) %>%
      summarise(across(all_of(to_plot_cols), sum)) %>%
      pivot_longer(cols = -all_of(agg_cols)) %>%
      mutate(
        type = name %>% str_remove_all("(_TP|_TN|_FP|_FN)"),
        measure = name %>% str_extract("(_TP|_TN|_FP|_FN)") %>% str_remove_all("_")
      ) %>%
      select(-name) %>%
      pivot_wider(id_cols = c(type, all_of(agg_cols)),
                  names_from = measure, values_from = value)
    
    if (to_plot == "sensitivity") {
      D = D %>% mutate(value = TP / (TP + FN))
    } else if (to_plot == "specificity") {
      D = D %>% mutate(value = TN / (TN + FP))
    }
  } else if (aggregate_prev & to_plot == "root mean square error") {
    D$AggPrev = D[[which_agg_prev]]
    agg_cols = c(facetting, x, "AggPrev")
    D = D %>%
      mutate(across(all_of(to_plot_cols), ~ (.x - AggPrev) ^ 2, .names = "{col}_SqErr")) %>%
      group_by(across(all_of(agg_cols))) %>%
      summarise(across(all_of(to_plot_cols), ~ sqrt(mean(.x, na.rm = T)), .names = "{col}")) %>%
      pivot_longer(cols = -all_of(agg_cols)) %>%
      mutate(type = name %>% str_remove_all("_ReproRate")) %>%
      select(-name)
    colnames(D)[colnames(D) == "AggPrev"] = which_agg_prev
  } else if (!aggregate_prev) {
    D = D %>% pivot_longer(cols = all_of(to_plot_cols), names_to = "type")
  }

  D
}

general_rep_plot = function (D, x, to_plot, facetting, types, to_include = NULL, show_points = F, aggregate_prev = F, which_agg_prev = "Prev_Lit", xlabel = NULL) {
  # browser()
  # Prepare tidy dataset for ggplot
  D = filtered_data(D, x, to_plot, facetting, types, to_include, aggregate_prev, which_agg_prev)
  D$xvar = D[[x]]
  
  # Build the plot
  if (length(facetting) == 1) {
    facet_formula = as.formula(paste("~", facetting[1]))
  } else if (length(facetting) == 2) {
    facet_formula = as.formula(paste(facetting[1], "~", facetting[2]))  
  }
  
  if (is.null(xlabel)) {
    xlabel = x
  }
  
  p = ggplot(D) +
    aes(x = xvar, y = value, group = type, color = type) +
    geom_line(stat = "summary", size = 0.65) +
    labs(
      color = "",
      x = xlabel,
      y = to_plot %>% str_to_sentence(),
      title = to_plot %>% str_to_sentence()
    )
  
  if (length(facetting) == 1) {
    p = p + facet_wrap(facet_formula)
  } else if (length(facetting) == 2) {
    p = p + facet_grid(facet_formula)
  }
  
  if (max(D$value <= 1, na.rm = T))
    p = p + scale_y_continuous(limits = c(0,1))
  else
    p = p + scale_y_continuous(limits = c(0,NA))
  
  if (is.numeric(D$xvar)) {
    if (max(D$xvar <= 1, na.rm = T))
      p = p + scale_x_continuous(limits = c(0,1))
    else
      p = p + scale_x_continuous(limits = c(0,NA))
  }
  
  if (show_points)
    p = p + geom_point(size = 1.5)
  
  if (to_plot == "reproducibility") 
    p = p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray")
  
  p
}

### Tracking plot
# X is prevalence, Y is reproducibility, sensitivity or specificity
# Lines represent the different measures of reproducibility (you can specify a subset)
# Facetted in grid by two of the parameters

tracking_plot = function (D, to_plot, prev, facetting, types, to_include = NULL, show_points = F, xlabel = NULL) {
  
  if (prev == "sample") x = "Prev_Sample"
  else if (prev == "literature") x = "Prev_Lit"

  general_rep_plot(D, x, to_plot, facetting, types, enquo(to_include), show_points, xlabel = xlabel)
}

### Aggregate tracking plot
# X is one of the parameters (except prevalence)
# Y is sensitivity or specificity
# Lines represent the different measures of reproducibility (you can specify a subset)
# Facetted in grid by two of the parameters (except prevalence)
# Data is aggregated for all prevalences available

aggregate_plot = function (D, x, to_plot, facetting, types, to_include = NULL, show_points = F, xlabel = NULL) {
  
  general_rep_plot(D, x, to_plot, facetting, types, enquo(to_include), show_points, aggregate_prev = T, which_agg_prev = "Prev_Lit", xlabel = xlabel)
 
}

### RMSE plot
# X is one of the parameters
# Y is RMSE
# Bar or lollipop plot
# Facetted in grid by two of the parameters

rmse_plot = function (D, x, prev, facetting, types, to_include = NULL, show_points = F, xlabel = NULL) {
  
  general_rep_plot(D, x, "root mean square error", facetting, types, enquo(to_include), show_points, aggregate_prev = T, which_agg_prev = prev, xlabel = xlabel)
  
}

plot_specification_curve = function (repro_rate_to_plot, spec_parameters, repro_rate_to_plot2 = NULL) {
  spec_param_values = map(spec_parameters, ~unique(DF[[.x]])) %>% unlist()
  
  DF_WIDE = bind_cols(
    map(1:length(spec_parameters), function (i) {
      par_col = spec_parameters[i]
      values = unique(DF[[par_col]])
      bind_cols(
        map(values, function (x) {
          df = data.frame(value = DF[[par_col]] == x)
          colnames(df) = x
          df
        })
      )
    })
  )
  
  DF = bind_cols(list(DF, DF_WIDE))
  # browser()
  DF$repro_rate = DF[[repro_rate_to_plot]]
  if (!is.null(repro_rate_to_plot2)) {
    DF$repro_rate2 = DF[[repro_rate_to_plot2]]
    DF_REPRO = DF %>% select(index, repro_rate, repro_rate2) %>%
      pivot_longer(cols = -index)
    plot_repro = ggplot(DF_REPRO) +
      aes(x = reorder(index, value), y = value, color = name) +
      theme(legend.position = "top")
  } else {
    plot_repro = ggplot(DF) +
      aes(x = reorder(index, repro_rate), y = repro_rate)
  }
  
  plot_repro = plot_repro +
    geom_point(size = 0.5) +
    labs(x = "", y = "Reproducibility\nRate") +
    ylim(c(0,1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), plot.title = element_blank())
  
  if (!is.null(repro_rate_to_plot2)) {
    cols_to_keep = c("index","repro_rate","repro_rate2", spec_param_values)
    SPECIFICATION = DF %>%
      select(all_of(cols_to_keep)) %>%
      pivot_longer(cols = -c(index, repro_rate, repro_rate2), names_to = "param")
  } else {
    cols_to_keep = c("index","repro_rate", spec_param_values)
    SPECIFICATION = DF %>%
      select(all_of(cols_to_keep)) %>%
      pivot_longer(cols = -c(index, repro_rate), names_to = "param")
  }
  
  plot_spec2 = ggplot(SPECIFICATION) +
    aes(x = reorder(index, repro_rate),
        y = param, fill = param, alpha = as.numeric(value)) +
    geom_tile(height = 0.4) +
    scale_alpha(range = c(0,1)) +
    # scale_fill_manual(breaks = spec_param_values, values = spec_colors) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "none", panel.grid = element_blank(),
          axis.text.y = element_text(size = 11), plot.title = element_blank())
  
  plot_grid(plot_repro, plot_spec2, align = "v", axis = "l", ncol = 1, nrow = 2, rel_heights = c(1,2))
}