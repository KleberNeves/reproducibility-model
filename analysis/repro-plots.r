setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")

# Actually getting the data
data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results/0720_2nd_Part_Run/"
repdata = get.figure.data.rep(data_dir)
repdata.backup = repdata

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName")
num_cols = !(colnames(repdata) %in% nonnum_cols)
repdata[,num_cols] = sapply(repdata[,num_cols], as.numeric)

repdata = repdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0("Power = ", 100 * typical.power, "%"),
         bias.label = paste0("Bias = ", 100 * bias.level, "%"))

repdata$scenarioName = car::recode(repdata$scenarioName, "'Peaks SD 0.5' = 'Overlapping Peaks'; 'Peaks SD 0.2' = 'Two Peaks'")


#########  Drawing reproducibility plots  ###########

figtheme = theme_linedraw()

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
