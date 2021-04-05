setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")
source("global.r")

# Load the data
data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results2/0421_Part2_All_Figures_Small2/"

# Prepare the data
repdata = get.figure.data.rep(data_dir)
repdata.backup = repdata

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

### Tracking plot
# X is prevalence, Y is reproducibility, sensitivity or specificity
# Lines represent the different measures of reproducibility (you can specify a subset)
# Facetted in grid by two of the parameters

tracking_plot = function (D, to_plot, prev, facetting, types, to_include = NULL) {
  
  # Filter the data if a filter is given
  if (!rlang::quo_is_null(enquo(to_include))) {
    D = D %>% filter(rlang::eval_tidy(enquo(to_include), D))
  }
  
  # Filter the types to plot
  if (to_plot == "reproducibility") { to_plot = "ReproRate" }
  else if (to_plot == "sensitivity") { to_plot = "SENS" }
  else if (to_plot == "specificity") { to_plot = "SPEC" }
  
  if (prev == "sample") prev = "Prev_Sample"
  else if (prev == "literature") prev = "Prev_Lit"
  
  to_plot = paste0(types, "_", to_plot)
  
  keep = c(to_plot, facetting, prev)
  D = D %>% select(all_of(keep))
  # browser()
  # Prepare tidy dataset for ggplot
  D = D %>%
    pivot_longer(cols = -all_of(c(prev, facetting))) %>%
    mutate(type = name %>% str_remove(paste0("_", to_plot)))
  
  D$prev = D[[prev]]
  # browser()
  # Build the plot
  p = ggplot(D) +
    aes(x = prev, y = value, group = type, color = type) +
    geom_line(stat = "summary") +
    geom_point() +
    facet_grid(as.formula(paste(facetting[1], "~", facetting[2])))
  
  p
}

types_to_plot = global_rep_types[c(1,2,9,12)]
tracking_plot(repdata, "reproducibility", prev = "literature", facetting = c("interlab.label", "power.label"), types = types_to_plot, to_include = (N == 10))



### RMSE plot
# X is either bias or power or interlab variation
# Y is reproducibility, sensitivity or specificity, averaged across prevalences
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

plot.rep.type = function (D, dist_shape, all_or_20, prev_all, which_measure, exclude_types = c()) {
  DS = D %>% filter(scenarioName == dist_shape &
                      Sample == all_or_20 &
                      Measure == which_measure &
                      !(Type %in% exclude_types))
  
  p = ggplot(DS)
  
  if (prev_all) {
    p = p + aes(x = `Prev_Sample_All`)
    xlabel = 'Prevalence of "True" Effects (literature)'
    prevlabel = "Prevalence (literature)"
  } else {
    p = p + aes(x = `Prev_Sample_20`)
    xlabel = 'Prevalence of "True" Effects (in the sample)'
    prevlabel = 'Prevalence (in the sample)'
  }
  
  tt = paste0(dist_shape, " - ", all_or_20, " exps, ", which_measure, " x ", prevlabel)
  
  p = p + aes(y = value, group = Type, color = as.factor(Type)) +
    geom_point(position = position_jitter(width = 0.0015)) +
    geom_line(stat = "summary", size = 0.6) +
    ylim(c(0,1)) + xlim(c(-0.05,1.05)) +
    facet_grid(bias.label ~ power.label) +
    scale_color_manual(values = c("brown", "springgreen4", "dodgerblue4", "dodgerblue2", "paleturquoise4")) +
    labs(x = xlabel, y = which_measure,
         title = tt, color = "Definition") +
    figtheme
  
  if (which_measure == "Reproducibility") {
    p = p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")
  }
  
  ggsave(plot = p,
         filename = paste0(tt, ".png"),
         width = 12, height = 7, dpi = 150)
  
  p
}

plot.error = function (D, dist_shape, all_or_20, prev_all, which_measure, which_compare) {
  DS = D %>% filter(scenarioName == dist_shape &
                      Sample == all_or_20 &
                      Measure == which_measure &
                      Compare == which_compare)
  
  p = ggplot(DS)
  
  if (prev_all) {
    p = p + aes(x = `Prev_Sample_All`)
    xlabel = 'Prevalence of "True" Effects (literature)'
    prevlabel = "Prevalence (literature)"
  } else {
    p = p + aes(x = `Prev_Sample_20`)
    xlabel = 'Prevalence of "True" Effects (in the sample)'
    prevlabel = 'Prevalence (in the sample)'
  }
  
  tt = paste0(dist_shape, " - ", all_or_20, " exps, ", which_measure, " (x ", which_compare,") x ", prevlabel)
  
  p = p + aes(y = value, group = bias.label, color = as.factor(bias.label)) +
    geom_point(position = position_jitter(width = 0.0015)) +
    geom_line(stat = "summary", size = 0.6) +
    xlim(c(-0.05,1.05)) +
    facet_wrap(~power.label) +
    scale_color_manual(values = c("blue","brown","green")) +
    labs(x = xlabel, y = paste0(which_measure, "(compared to ", which_compare, ")"),
         title = tt, color = "Bias") +
    figtheme
  
  if (which_measure != "Exaggeration") {
    p = p + ylim(c(0,1))
  }
  
  ggsave(plot = p,
         filename = paste0(tt, ".png"),
         width = 12, height = 4, dpi = 150)
  
  p
}

# Plots the reproducibility, sensitivity or specificity of the many definitions of "success",
# facetting for power and interlab variation (instead of bias)
plot.rep.type2 = function (D, dist_shape, all_or_20, prev_all, which_measure, exclude_types = c()) {
  DS = D %>% filter(scenarioName == dist_shape &
                      Sample == all_or_20 &
                      Measure == which_measure &
                      !(Type %in% exclude_types))
  
  p = ggplot(DS)
  
  if (prev_all) {
    p = p + aes(x = `Prev_SampleBias_All`)
    xlabel = 'Prevalence of Unbiased Effects (literature)'
    prevlabel = "Prevalence (literature)"
  } else {
    p = p + aes(x = `Prev_SampleBias_20`)
    xlabel = 'Prevalence of Unbiased Effects (in the sample)'
    prevlabel = 'Prevalence (in the sample)'
  }
  
  tt = paste0(dist_shape, " - ", all_or_20, " exps, ", which_measure, " x ", prevlabel)
  
  p = p + aes(y = value, group = Type, color = as.factor(Type)) +
    geom_point() +
    geom_line(stat = "summary", size = 0.6) +
    ylim(c(0,1)) + xlim(c(-0.05,1.05)) +
    facet_grid(interlab.var.label ~ power.label) +
    scale_color_manual(values = c("brown", "springgreen4", "dodgerblue4", "dodgerblue2", "paleturquoise4")) +
    labs(x = xlabel, y = which_measure,
         title = tt, color = "Definition") +
    figtheme
  
  if (which_measure == "Reproducibility") {
    p = p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey")
  }
  
  ggsave(plot = p,
         filename = paste0(tt, ".png"),
         width = 12, height = 7, dpi = 150)
  
  p
}

### REPRODUCIBILITY RATE PLOTS

D = repdata %>% select(scenarioName, Prev_Sample_20, Prev_Sample_All, bias.label, power.label, RepSet, 21:50)
D = D %>% melt(id.vars = c("scenarioName", "Prev_Sample_20", "Prev_Sample_All", "bias.label", "power.label", "RepSet"))
D$Sample = ifelse(str_detect(D$variable, "_All"), "All", "20")
D$Measure = str_remove_all(str_extract(D$variable, "_.+?_"),"_")
D$Measure = car::recode(D$Measure, "'ReproRate' = 'Reproducibility'")
D$Type = str_remove_all(str_extract(D$variable, ".+?_"),"_")

plot.rep.type(D, "Single Normal", "All", T, "Reproducibility")
plot.rep.type(D, "Single Normal", "20", T, "Reproducibility")
plot.rep.type(D, "Single Normal", "20", F, "Reproducibility")
plot.rep.type(D, "Single Normal", "All", F, "Specificity")
plot.rep.type(D, "Single Normal", "20", F, "Specificity")
plot.rep.type(D, "Single Normal", "All", F, "Sensitivity")
plot.rep.type(D, "Single Normal", "20", F, "Sensitivity")

plot.rep.type(D, "Dichotomous", "All", T, "Reproducibility")
plot.rep.type(D, "Dichotomous", "20", T, "Reproducibility")
plot.rep.type(D, "Dichotomous", "20", F, "Reproducibility")
plot.rep.type(D, "Dichotomous", "All", F, "Specificity")
plot.rep.type(D, "Dichotomous", "20", F, "Specificity")
plot.rep.type(D, "Dichotomous", "All", F, "Sensitivity")
plot.rep.type(D, "Dichotomous", "20", F, "Sensitivity")


### ERROR PLOTS

D = repdata %>% select(scenarioName, Prev_Sample_20, Prev_Sample_All, bias.label, power.label, RepSet, 51:58)
D = D %>% melt(id.vars = c("scenarioName", "Prev_Sample_20", "Prev_Sample_All", "bias.label", "power.label", "RepSet"))
D$Sample = ifelse(str_detect(D$variable, "_All"), "All", "20")
D$Measure = str_remove_all(str_extract(D$variable, ".+?_"),"_")
D$Compare = str_remove_all(str_extract(D$variable, "x_.+?_"),"[_x]")

plot.error(D, "Single Normal", "All", F, "Exaggeration", "Original")
plot.error(D, "Single Normal", "20", F, "Exaggeration", "Real")
plot.error(D, "Single Normal", "All", F, "Signal", "Original")
plot.error(D, "Single Normal", "20", F, "Signal", "Real")

plot.error(D, "Dichotomous", "All", F, "Exaggeration", "Original")
plot.error(D, "Dichotomous", "20", F, "Exaggeration", "Real")
plot.error(D, "Dichotomous", "All", F, "Signal", "Original")
plot.error(D, "Dichotomous", "20", F, "Signal", "Real")


####### REPRODUCIBILITY PLOTS BIAS-BASED

D = repdata %>% select(scenarioName, Prev_SampleBias_20, Prev_SampleBias_All, interlab.var.label, bias.label, power.label, RepSet, 21:70)
D = D %>% melt(id.vars = c("scenarioName", "Prev_SampleBias_20", "Prev_SampleBias_All", "interlab.var.label", "bias.label", "power.label", "RepSet"))
D$Sample = ifelse(str_detect(D$variable, "_All"), "All", "20")
D$Measure = str_remove_all(str_extract(D$variable, "_.+?_"),"_")
D$Measure = car::recode(D$Measure, "'ReproRate' = 'Reproducibility'")
D$Type = str_remove_all(str_extract(D$variable, ".+?_"),"_")

plot.rep.type2(D, "Two Peaks (SD = 0.1)", "All", T, "Reproducibility")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "20", T, "Reproducibility")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "20", F, "Reproducibility")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "All", F, "Specificity")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "All", F, "SpecificityBias")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "All", F, "Sensitivity")
plot.rep.type2(D, "Two Peaks (SD = 0.1)", "All", F, "SensitivityBias")