setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")

# Actually getting the data
data_dir = "../../../Results/0720_2nd_Part_Test/"
figdata = get.figure.data.eval(data_dir)

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName")
num_cols = !(colnames(figdata) %in% nonnum_cols)
figdata[,num_cols] = sapply(figdata[,num_cols], as.numeric)

figdata = figdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0("Power = ", 100 * typical.power, "%"),
         bias.label = paste0("Bias = ", 100 * bias.level, "%"))

figdata$scenarioName = car::recode(figdata$scenarioName, "'Peaks SD 0.5' = 'Overlapping Peaks'; 'Peaks SD 0.2' = 'Two Peaks'")


#########  Drawing plots  ###########

figtheme = theme_linedraw()

# Figures similar to Ioannidis (2005), lines are bias, using odds
D = figdata %>% filter(true_odds < 1.5)

ps = dlply(D, "scenarioName", function (D) {
  p = ggplot(D) +
    aes(x = true_odds, group = bias.level, color = as.factor(bias.level)) +
    geom_line(size = 0.7) +
    facet_wrap(~power.label) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("blue","brown","green")) +
    labs(x = 'Odds Ratio for "True" Effects',
         title = D$scenarioName[1], color = "Bias") +
    figtheme
  
  pp = p + aes(y = `Positive Predictive Value (Correct Signal)_Rate`) +
    labs(y = "Positive Predictive Value") + ylim(c(0,1))
  
  ps = p + aes(y = `Signal Error_Rate`) +
     labs(y = "Signal Error Rate") + ylim(c(0,1))
  
  pm = p + aes(y = `Exaggeration Factor_Median`) +
    aes(ymin = `Exaggeration Factor_Q25`,
        ymax = `Exaggeration Factor_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio")
  
  ps2 = p + aes(y = `Signal Error (Effects > Min)_Rate`) +
    labs(y = "Signal Error Rate (for true positives only)") + ylim(c(0,1))
  
  pm2 = p + aes(y = `Exaggeration Factor (Effects > Min)_Median`) +
    aes(ymin = `Exaggeration Factor (Effects > Min)_Q25`,
        ymax = `Exaggeration Factor (Effects > Min)_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio (for true positives only)")
  
  ggsave(plot = pp,
         filename = paste0("By Power, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- PPV.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Power, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Power, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps2,
         filename = paste0("By Power, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm2,
         filename = paste0("By Power, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  list(pp, ps, pm, ps2, pm2)
})

# Figures similar to Ioannidis (2005), lines are bias, using prevalence
D = figdata

ps = dlply(D, "scenarioName", function (D) {
  p = ggplot(D) +
    aes(x = `% Above minimum`, group = bias.level, color = as.factor(bias.level)) +
    geom_line(size = 0.7) +
    facet_wrap(~power.label) +
    scale_color_manual(values = c("blue","brown","green")) +
    labs(x = 'Prevalence of "True" Effects',
         title = D$scenarioName[1], color = "Bias") +
    figtheme
  
  pp = p + aes(y = `Positive Predictive Value (Correct Signal)_Rate`) +
    labs(y = "Positive Predictive Value") + ylim(c(0,1))
  
  ps = p + aes(y = `Signal Error_Rate`) +
    labs(y = "Signal Error Rate") + ylim(c(0,1))
  
  pm = p + aes(y = `Exaggeration Factor_Median`) +
    aes(ymin = `Exaggeration Factor_Q25`,
        ymax = `Exaggeration Factor_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio")
  
  ps2 = p + aes(y = `Signal Error (Effects > Min)_Rate`) +
    labs(y = "Signal Error Rate (for true positives only)") + ylim(c(0,1))
  
  pm2 = p + aes(y = `Exaggeration Factor (Effects > Min)_Median`) +
    aes(ymin = `Exaggeration Factor (Effects > Min)_Q25`,
        ymax = `Exaggeration Factor (Effects > Min)_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio (for true positives only)")
  
  ggsave(plot = pp,
         filename = paste0("By Power, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- PPV.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Power, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Power, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps2,
         filename = paste0("By Power, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm2,
         filename = paste0("By Power, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  list(pp, ps, pm, ps2, pm2)
})


# Figures similar to Ioannidis (2005), but lines are power, using odds
D = figdata %>% filter(true_odds < 1.5)

ps = dlply(D, "scenarioName", function (D) {
  p = ggplot(D) +
    aes(x = true_odds, group = power.label, color = as.factor(power.label)) +
    geom_line(size = 0.7) +
    facet_wrap(~bias.label) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("blue","brown","green")) +
    labs(x = 'Odds Ratio for "True" Effects',
         title = D$scenarioName[1], color = "Power") +
    figtheme
  
  pp = p + aes(y = `Positive Predictive Value (Correct Signal)_Rate`) +
    labs(y = "Positive Predictive Value") + ylim(c(0,1))
  
  ps = p + aes(y = `Signal Error_Rate`) +
    labs(y = "Signal Error Rate") + ylim(c(0,1))
  
  pm = p + aes(y = `Exaggeration Factor_Median`) +
    aes(ymin = `Exaggeration Factor_Q25`,
        ymax = `Exaggeration Factor_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio")
  
  ps2 = p + aes(y = `Signal Error (Effects > Min)_Rate`) +
    labs(y = "Signal Error Rate (for true positives only)") + ylim(c(0,1))
  
  pm2 = p + aes(y = `Exaggeration Factor (Effects > Min)_Median`) +
    aes(ymin = `Exaggeration Factor (Effects > Min)_Q25`,
        ymax = `Exaggeration Factor (Effects > Min)_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio (for true positives only)")
  
  ggsave(plot = pp,
         filename = paste0("By Bias, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- PPV.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Bias, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Bias, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Bias, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Bias, X Axis Odds - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  list(pp, ps, pm, ps2, pm2)
})

# Figures similar to Ioannidis (2005), but lines are power, using prevalence
D = figdata

ps = dlply(D, "scenarioName", function (D) {
  p = ggplot(D) +
    aes(x = `% Above minimum`, group = power.label, color = as.factor(power.label)) +
    geom_line(size = 0.7) +
    facet_wrap(~bias.label) +
    scale_color_manual(values = c("blue","brown","green")) +
    labs(x = 'Prevalence of "True" Effects',
         title = D$scenarioName[1], color = "Power") +
    figtheme
  
  pp = p + aes(y = `Positive Predictive Value (Correct Signal)_Rate`) +
    labs(y = "Positive Predictive Value") + ylim(c(0,1))
  
  ps = p + aes(y = `Signal Error_Rate`) +
    labs(y = "Signal Error Rate") + ylim(c(0,1))
  
  pm = p + aes(y = `Exaggeration Factor_Median`) +
    aes(ymin = `Exaggeration Factor_Q25`,
        ymax = `Exaggeration Factor_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio")
  
  ps2 = p + aes(y = `Signal Error (Effects > Min)_Rate`) +
    labs(y = "Signal Error Rate (for true positives only)") + ylim(c(0,1))
  
  pm2 = p + aes(y = `Exaggeration Factor (Effects > Min)_Median`) +
    aes(ymin = `Exaggeration Factor (Effects > Min)_Q25`,
        ymax = `Exaggeration Factor (Effects > Min)_Q75`) + geom_errorbar() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(y = "Exaggeration Ratio (for true positives only)")
  
  ggsave(plot = pp,
         filename = paste0("By Bias, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- PPV.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Bias, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Bias, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M.png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = ps,
         filename = paste0("By Bias, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type S (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  ggsave(plot = pm,
         filename = paste0("By Bias, X Axis Prev - ",
                           str_replace_all(D$scenarioName[1], "[.]", "_"), "- Type M (for true positives only).png"),
         width = 12, height = 4, dpi = 150)
  
  list(pp, ps, pm, ps2, pm2)
})

