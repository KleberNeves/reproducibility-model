setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")
source("global.r")

# Load the data
data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results2/0421_Part2_Sweep1/"

# Prepare the data

# To process the data by parts
zips = paste(data_dir, list.files(data_dir, ".zip$"), sep ="")
total = length(zips)
i = 0
process_df = tibble(i = 1:total, processed = F)
repdata_list = vector(mode = "list", length = total)

walk(zips, function (z) {
  i <<- i + 1
  print(paste0(i, "/", total))
  get.data.from.zip.rep(z)
})


# To process the data at once
# repdata = get.figure.data.rep(data_dir)


# repdata.backup = repdata
# save(repdata, file = paste0("repdata.RData"))
# To load previously processed data
load("repdata.RData")

repdata = repdata %>%
  mutate(across(where(~ !anyNA(as.numeric(.x))), as.numeric))

repdata = repdata %>%
  mutate(
    power.label = paste0("Power = ", 100 * typical.power, "%"),
    bias.label = paste0("Bias = ", 100 * bias.level, "%"),
    interlab.label = paste0("Interlab Var = ", round(100 * interlab.var / (interlab.var + 1)), "%"),
  )

repdata$scenarioName = car::recode(repdata$scenarioName, "'Peaks SD 0.5' = 'Overlapping Peaks'; 'Peaks SD 0.2' = 'Two Peaks'")


##### Plot functions #####

figtheme = theme_linedraw() +
  theme(
    legend.position = "bottom"
  )

theme_set(figtheme)

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
      pivot_wider(id_cols = c(type, all_of(agg_cols)), names_from = measure, values_from = value)
    
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
  }

  D
}

general_rep_plot = function (D, x, to_plot, facetting, types, to_include = NULL, show_points = F, aggregate_prev = F, which_agg_prev = "Prev_Lit", xlabel = NULL) {
  
  # Prepare tidy dataset for ggplot
  D = filtered_data(D, x, to_plot, facetting, types, to_include, aggregate_prev, which_agg_prev)
  D$xvar = D[[x]]
  
  # Build the plot
  facet_formula = as.formula(paste(facetting[1], "~", facetting[2]))
  if (is.null(xlabel)) {
    xlabel = x
  }
  
  p = ggplot(D) +
    aes(x = xvar, y = value, group = type, color = type) +
    geom_line(stat = "summary", size = 0.65) +
    facet_grid(facet_formula) +
    labs(
      color = "",
      x = xlabel,
      y = to_plot %>% str_to_sentence(),
      title = to_plot %>% str_to_sentence()
    )
  
  if (max(D$value <= 1))
    p = p + scale_y_continuous(limits = c(0,1))
  else
    p = p + scale_y_continuous(limits = c(0,NA))
  
  if (is.numeric(D$xvar)) {
    if (max(D$value <= 1))
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

types_to_plot = global_rep_types[c(1,2,9,12)]
rmse_plot(repdata, x = "interlab.label", prev = "Prev_Lit", facetting = c("bias.label", "power.label"), types = types_to_plot, to_include = (N == 10))

##### Test plots #####

types_to_plot = global_rep_types[c(1,2,9,12)]
tracking_plot(repdata, "reproducibility", prev = "literature", facetting = c("interlab.label", "power.label"), types = types_to_plot, to_include = (N == 10))

types_to_plot = global_rep_types[c(1,2,9,12)]
aggregate_plot(repdata, x = "interlab.label", "sensitivity", facetting = c("bias.label", "power.label"), types = types_to_plot, to_include = (N == 10))
