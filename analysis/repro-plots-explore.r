setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")
source("global.r")

# Load the data
data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results2/0421_Part2_All_Figures_Small2/"

# Prepare the data
repdata = get.figure.data.rep(data_dir)
# repdata.backup = repdata
save(repdata, file = paste0("repdata.RData"))
# load("repdata.RData")

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

filtered_data = function (D, x, to_plot, facetting, types, to_include, aggregate_prev = F) {
  # Filter the data if a filter is given
  if (!rlang::quo_is_null(to_include)) {
    D = D %>% filter(rlang::eval_tidy(to_include, D))
  }
  
  # Filter the types to plot
  if (to_plot == "reproducibility") { to_plot_suffix = "ReproRate" }
  else if (to_plot == "sensitivity") { to_plot_suffix = "SENS" }
  else if (to_plot == "specificity") { to_plot_suffix = "SPEC" }
  
  to_plot_cols = paste0(types, "_", to_plot_suffix)
  
  keep = c(to_plot_cols, facetting, x)
  if (aggregate_prev) keep = c(keep, "Prev_Lit")
  
  D = D %>% select(all_of(keep))
  browser()
  if (aggregate_prev) {
    agg_cols = c(facetting, x, "Prev_Lit")
    D = D %>% group_by(across(all_of(agg_cols))) %>%
      summarise(across(all_of(to_plot_cols), sum))
  }
  
  D = D %>%
    pivot_longer(cols = -all_of(c(x, facetting))) %>%
    mutate(type = name %>% str_remove(paste0("_", to_plot_suffix)))
}

general_rep_plot = function (D, x, to_plot, facetting, types, to_include = NULL, show_points = F, aggregate_prev = F) {
  
  # Prepare tidy dataset for ggplot
  D = filtered_data(D, x, to_plot, facetting, types, to_include, aggregate_prev)
  D$xvar = D[[x]]
  
  # Build the plot
  facet_formula = as.formula(paste(facetting[1], "~", facetting[2]))
  
  p = ggplot(D) +
    aes(x = xvar, y = value, group = type, color = type) +
    geom_line(stat = "summary", size = 0.65) +
    scale_y_continuous(limits = c(0,1)) +
    facet_grid(facet_formula) +
    labs(
      color = "",
      x = x,
      y = to_plot %>% str_to_title()
    )
  
  if (is.numeric(D$xvar))
    p = p + scale_x_continuous(limits = c(0,1))
  
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

tracking_plot = function (D, to_plot, prev, facetting, types, to_include = NULL, show_points = F) {
  
  if (prev == "sample") x = "Prev_Sample"
  else if (prev == "literature") x = "Prev_Lit"

  general_rep_plot(D, x, to_plot, facetting, types, enquo(to_include), show_points)
}

### Aggregate tracking plot
# X is one of the parameters (except prevalence)
# Y is sensitivity or specificity
# Lines represent the different measures of reproducibility (you can specify a subset)
# Facetted in grid by two of the parameters (except prevalence)
# Data is aggregated for all prevalences available

aggregate_plot = function (D, x, to_plot, facetting, types, to_include = NULL, show_points = F) {
  
  general_rep_plot(D, x, to_plot, facetting, types, enquo(to_include), show_points, aggregate_prev = T)
 
}

aggregate_plot(repdata, x = "interlab.label", "sensitivity", facetting = c("bias.label", "power.label"), types = types_to_plot, to_include = (N == 10))

### RMSE plot
# X is one of the parameters
# Y is RMSE
# Bar or lollipop plot
# Facetted in grid by two of the parameters



### ROC Curve plot
# X is FP / (FP + TN)
# Y is TP / (TP + FN)
# Points
# Color represent either bias, power or interlab variation
# Facetted in grid by two of the parameters

### Extra settings
# For facetting or coloring or X axis, parameters can be:
#   bias, power, interlab variation, repro.power, etc
# Additionally, you can pass a filter to the functions, to be applied before
# All these setting are registered in the caption of the figure and in the filename

##### Test plots #####

types_to_plot = global_rep_types[c(1,2,9,12)]
tracking_plot(repdata, "reproducibility", prev = "literature", facetting = c("interlab.label", "power.label"), types = types_to_plot, to_include = (N == 10))

aggregate_plot(repdata, "sensitivity", prev = "literature", x = "interlab.label", facetting = c("bias.label", "power.label"), types = types_to_plot, to_include = (N == 10))