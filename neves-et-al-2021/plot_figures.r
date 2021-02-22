##### Configuration #####

path_data_main_sims = "./Data/MainFigures/"
path_data_only_if_large_sims = "./Data/PublishOnlyIfLarge/"
path_data_power_calc_sims = "./Data/PowerCalcVariations/"
path_data_table_sims = "./Data/Table/"


# Figures will be saved in the working directory
path_output = "."

# If running on RStudio, use the line below to set working directory to this file's path
path_output = setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### Script #####

source("helper.r")
source("model.r")

library(cowplot)
library(purrr)
library(glue)

dir.create(paste0(path_output, "/Main"))
dir.create(paste0(path_output, "/Publish only if large"))
dir.create(paste0(path_output, "/Power calc variations"))

##### Function definition #####

# Functions for plotting

plot_figure = function(ps, variable, fig_title, models) {
  plots = map(ps, function(x) { x[[variable]] })
  plots = plots[models]
  fig = plot_grid(plotlist = plots, ncol = 1, labels = "AUTO", label_size = 14)
  ggsave2(glue("{fig_title}.png"), fig, width = 9, height = 3 * length(plots), dpi = 150)
}

plot.reference.distribution = function (a.input, fname) {
  uldist = data.frame(x = sample.from.dist(a.input$sdA, a.input$weightB,
                                           a.input$meanB, a.input$sdB, k = 5000000
  ))
  
  spread = min(max(a.input$sdA, a.input$sdB, 0.5), 1) * 4
  
  p = ggplot(uldist, aes(x = x)) +
    geom_vline(xintercept = a.input$m,
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = -a.input$m,
               linetype = "dashed", color = "blue") +
    labs(title = fname, y = "", x = "") + xlim(c(-spread,spread)) +
    theme_minimal_grid() +
    theme(
      axis.text.y = element_blank(),
      plot.title = element_text(size = 10)
    )
  
  if (a.input$sdB == a.input$sdA & a.input$sdA == 0) {
    p = p + geom_histogram()
  } else {
    p = p + geom_density()
  }
  
  p
}

make_plots = function (D) {
  ps = dlply(D, "scenarioName", function (D) {
    p = ggplot(D) +
      aes(x = `% Above minimum`, group = bias.level, color = as.factor(bias.level)) +
      geom_line(size = 0.7) +
      geom_point(size = 1) +
      facet_wrap(~power.label) +
      scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(0,1)) +
      scale_color_manual(values = c("blue","brown","green")) +
      labs(x = 'Prevalence of True Effects',
           title = D$scenarioName[1], color = "Bias") +
      figtheme
    
    pp = p + aes(y = `Positive Predictive Value Signal_Rate`) +
      labs(y = "Positive Predictive Value") + scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,1))
    
    ps = p + aes(y = `Signal Error_Rate`) +
      labs(y = "Signal Error Rate") + scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,0.5))
    
    pm = p + aes(y = `Exaggeration Factor_Median`) +
      aes(ymin = `Exaggeration Factor_Q25`,
          ymax = `Exaggeration Factor_Q75`) + geom_errorbar(width = 0) +
      scale_y_continuous(breaks = scales::pretty_breaks(), trans = "log10") +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Exaggeration Ratio")
    
    ps2 = p + aes(y = `Signal Error (Effects > Min)_Rate`) +
      labs(y = "Signal Error Rate (for true positives only)") + scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,0.5))
    
    pm2 = p + aes(y = `Exaggeration Factor (Effects > Min)_Median`) +
      aes(ymin = `Exaggeration Factor (Effects > Min)_Q25`,
          ymax = `Exaggeration Factor (Effects > Min)_Q75`) + geom_errorbar(width = 0) +
      scale_y_continuous(breaks = scales::pretty_breaks(), trans = "log10") +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Exaggeration Ratio (for true positives only)")
    
    list(PPV = pp, Signal = ps, Magnitude = pm, Signal_True = ps2, Magnitude_True = pm2)
  })
  
  ps
}

make_plots_power_calculation = function (D) {
  ps = dlply(D, "power.calculation", function (D) {
    p = ggplot(D) +
      aes(x = `% Above minimum`, group = bias.level, color = as.factor(bias.level)) +
      geom_line(size = 0.7) +
      geom_point(size = 1) +
      facet_wrap(~power.label) +
      scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(0,1)) +
      scale_color_manual(values = c("blue","brown","green")) +
      labs(x = 'Prevalence of True Effects',
           title = D$power.calculation[1], color = "Bias") +
      figtheme
    
    pp = p + aes(y = `Positive Predictive Value Signal_Rate`) +
      labs(y = "Positive Predictive Value") + scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,1))
    
    ps = p + aes(y = `Signal Error_Rate`) +
      labs(y = "Signal Error Rate") + scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,0.5))
    
    pm = p + aes(y = `Exaggeration Factor_Median`) +
      aes(ymin = `Exaggeration Factor_Q25`,
          ymax = `Exaggeration Factor_Q75`) + geom_errorbar(width = 0) +
      scale_y_continuous(breaks = c(2,4,6,8,12), limits = c(1, 13), trans = "log10") +
      # geom_hline(yintercept = 1, linetype = "dashed") +
      labs(y = "Exaggeration Ratio")
    
    list(PPV = pp, Signal = ps, Magnitude = pm)
  })
  
  ps
}

##### ggplot theme #####

figtheme = theme_minimal_grid() +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13),
    axis.title.x = element_text(size = 14, margin = margin(8,0,0,0)),
    axis.title.y = element_text(size = 12, margin = margin(0,8,0,0)),
    plot.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.spacing.x = unit(16,"points")
  )

##### Main simulations #####
data_dir = path_data_main_sims
figdata = get.figure.data.eval(data_dir, only_published = T)

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName", "power.calculation")
num_cols = !(colnames(figdata) %in% nonnum_cols)
figdata[,num_cols] = sapply(figdata[,num_cols], as.numeric)

figdata = figdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0("Power = ", 100 * typical.power, "%"),
         bias.label = paste0("Bias = ", 100 * bias.level, "%"))

figdata$scenarioName = car::recode(figdata$scenarioName, "'Peaks SD 0.3' = 'Overlapping Peaks'; 'Peaks SD 0.1' = 'Two Peaks'")

##### Plot main paper figures #####
D = figdata

ps = make_plots(D)

all_models = c("Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")
all_continuous_models = c("Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")
main_models = c("Two Peaks", "Single Normal")
other_models = c("Overlapping Peaks", "Two Normals")

plot_figure(ps, variable = "PPV",
            fig_title = paste0(path_output, "/Main/Figure - PPV"),
            models = all_models)

plot_figure(ps, variable = "Signal",
            fig_title = paste0(path_output, "/Main/Figure - Signal"),
            models = all_continuous_models)

plot_figure(ps, variable = "Magnitude",
            fig_title = paste0(path_output, "/Main/Figure - Magnitude"),
            models = all_continuous_models)

# Odds-axis figure to be compared to Ioannidis 2005

p = ggplot(D %>% filter(true_odds < 1.5, scenarioName == "Dichotomous")) +
  aes(x = true_odds, y = `Positive Predictive Value_Rate`,
      group = bias.level, color = as.factor(bias.level)) +
  geom_line(size = 0.7) +
  geom_point(size = 1) +
  facet_wrap(~power.label) +
  scale_x_continuous(breaks = scales::pretty_breaks(), limits = c(0,1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0,1)) +
  scale_color_manual(values = c("blue","brown","green")) +
  labs(x = 'Odds of True Effects', y = "Positive Predictive Value",
       title = D$scenarioName[1], color = "Bias") +
  figtheme

ggsave2(paste0(path_output, "/Main/Figure - Odds PPV.png"), p, width = 9, height = 3, dpi = 150)


# Reference distribution plots

pE = plot.reference.distribution(
  list(sdA = 1, weightB = 0, meanB = 0, sdB = 0, m = 0.5), "Single Normal")
pD = plot.reference.distribution(
  list(sdA = 0.1, weightB = 0.5, meanB = 0, sdB = 1, m = 0.5), "Two Normals")
pC = plot.reference.distribution(
  list(sdA = 0.3, weightB = 0.5, meanB = 1, sdB = 0.3, m = 0.5), "Overlapping Peaks")
pB = plot.reference.distribution(
  list(sdA = 0.1, weightB = 0.5, meanB = 1, sdB = 0.1, m = 0.5), "Two Peaks")
pA = plot.reference.distribution(
  list(sdA = 0, weightB = 0.5, meanB = 1, sdB = 0, m = 0.5), "Dichotomous")

top_row = plot_grid(pA, pB, pC, ncol = 3, nrow = 1, labels = c("A","B","C"), label_size = 14)
bottom_row = plot_grid(pD, pE, ncol = 2, nrow = 1, labels = c("D","E"), label_size = 14)

fig = plot_grid(top_row, bottom_row, ncol = 1, nrow = 2)

ggsave2(paste0(path_output, "/Main/Figure - Distributions.png"), fig, width = 8, height = 4, dpi = 150)

##### Table #####
data_dir = path_data_table_sims
figdata = get.figure.data.eval(data_dir)

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName")
num_cols = !(colnames(figdata) %in% nonnum_cols)
figdata[,num_cols] = sapply(figdata[,num_cols], as.numeric)

figdata = figdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0(100 * typical.power, "%"),
         bias.label = paste0(100 * bias.level, "%"))

# Build the table
D = figdata

D$UnderlyingDist = str_remove(str_extract(D$scenarioName, ".+;"),";")
D$Scenario = str_remove(str_extract(D$scenarioName, ";.+"),"; ")
D$Scenario = factor(D$Scenario, levels = c("Adequately powered RCT with little bias and 1:1 pre-study odds", "Confirmatory meta-analysis of good quality RCTs", "Meta-analysis of small inconclusive studies", "Underpowered, but well-performed phase I/II RCT", "Underpowered, poorly performed phase I/II RCT", "Adequately powered exploratory epidemiological study", "Underpowered exploratory epidemiological study", "Discovery-oriented exploratory research with massive testing", "Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"))

# PPV

TB1 = D %>% select(power.label, bias.label, `Positive Predictive Value Signal_Rate`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Overlapping Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")

TB1P = TB1

# Type S

TB1 = D %>% select(power.label, bias.label, `Signal Error_Rate`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Overlapping Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")

TB1S = TB1

# Type M

TB1 = D %>% select(power.label, bias.label, `Exaggeration Factor_Median`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Overlapping Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")

TB1M = TB1

# Table comparing percentage above minimum of interest with the prevalence

TB2 = D %>% select(`% Above minimum`, Scenario, UnderlyingDist)
TB2 = TB2 %>% melt(id.vars = c("Scenario", "UnderlyingDist"))

TB2 = TB2 %>% select(-variable) %>%
  dcast(Scenario ~ UnderlyingDist, value.var = "value")

TB2$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB2$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB2 = TB2 %>% select(Scenario, Odds, Prev, Dichotomous, `Two Peaks`, `Overlapping Peaks`, `Two Normals`, `Single Normal`)

colnames(TB2) = c("Scenario description", "Odds for “True” Effects", "Prevalence of “True” Effects", "Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")


# Writing the tables to Excel
library(openxlsx)

# Adds column for spacing
TB1P = TB1P %>% mutate(` ` = "") %>% select(1:5,11,6:10)
TB1S = TB1S %>% mutate(` ` = "") %>% select(1:5,11,6:10)
TB1M = TB1M %>% mutate(` ` = "") %>% select(1:5,11,6:10)
TB2 = TB2 %>% mutate(` ` = "") %>% select(1:3,9,4:8)

wb = createWorkbook()

addWorksheet(wb, "PPV")
sht = 1
writeData(wb, sheet = sht, x = TB1P)
setColWidths(wb, sht, cols = 1:11, widths = c(60,10,10,10,10,3,12,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:11, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:11, gridExpand = TRUE)

addWorksheet(wb, "Type S")
sht = 2
writeData(wb, sheet = sht, x = TB1S)
setColWidths(wb, sht, cols = 1:11, widths = c(60,10,10,10,10,3,12,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:11, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:11, gridExpand = TRUE)

addWorksheet(wb, "Type M")
sht = 3
writeData(wb, sheet = sht, x = TB1M)
setColWidths(wb, sht, cols = 1:11, widths = c(60,10,10,10,10,3,12,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:11, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:11, gridExpand = TRUE)

addWorksheet(wb, "Prevalence")
sht = 4
writeData(wb, sheet = sht, x = TB2)
setColWidths(wb, sht, cols = 1:9, widths = c(60,10,10,3,12,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:9, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:9, gridExpand = TRUE)

saveWorkbook(wb, paste0(path_output, "/Main/Table Comparison Ioannidis.xlsx"), overwrite = TRUE)

##### "Publish only if large" simulations #####
data_dir = path_data_only_if_large_sims
figdata = get.figure.data.eval(data_dir, only_published = T)

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName", "power.calculation")
num_cols = !(colnames(figdata) %in% nonnum_cols)
figdata[,num_cols] = sapply(figdata[,num_cols], as.numeric)

figdata = figdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0("Power = ", 100 * typical.power, "%"),
         bias.label = paste0("Bias = ", 100 * bias.level, "%"))

figdata$scenarioName = car::recode(figdata$scenarioName, "'Peaks SD 0.3' = 'Overlapping Peaks'; 'Peaks SD 0.1' = 'Two Peaks'")

##### Plot supplementary figures #####
D = figdata

ps = make_plots(D)

all_models = c("Dichotomous", "Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")
all_continuous_models = c("Two Peaks", "Overlapping Peaks", "Two Normals", "Single Normal")
main_models = c("Two Peaks", "Single Normal")
other_models = c("Overlapping Peaks", "Two Normals")

plot_figure(ps, variable = "PPV",
            fig_title = paste0(path_output, "/Publish only if large/Figure - PPV"),
            models = all_models)

plot_figure(ps, variable = "Signal",
            fig_title = paste0(path_output, "/Publish only if large/Figure - Signal"),
            models = all_continuous_models)

plot_figure(ps, variable = "Magnitude",
            fig_title = paste0(path_output, "/Publish only if large/Figure - Magnitude"),
            models = all_continuous_models)

##### Power calculation variations simulations #####
data_dir = path_data_power_calc_sims
figdata = get.figure.data.eval(data_dir, only_published = T)

# Convert to numeric
nonnum_cols = c("calc.repro", "how.sim.ends", "scenarioName", "power.calculation")
num_cols = !(colnames(figdata) %in% nonnum_cols)
figdata[,num_cols] = sapply(figdata[,num_cols], as.numeric)

figdata = figdata %>%
  mutate(param_odds = weightB / (1 - weightB),
         true_odds = `% Above minimum` / (1 - `% Above minimum`),
         power.label = paste0("Power = ", 100 * typical.power, "%"),
         bias.label = paste0("Bias = ", 100 * bias.level, "%"))

figdata$scenarioName = car::recode(figdata$scenarioName, "'Peaks SD 0.3' = 'Overlapping Peaks'; 'Peaks SD 0.1' = 'Two Peaks'")

figdata$power.calculation = str_to_sentence(figdata$power.calculation)

ps = make_plots_power_calculation(figdata)

all_models = unique(figdata$power.calculation)

plot_figure(ps, variable = "PPV",
            fig_title = paste0(path_output, "/Power calc variations/Figure - PPV"),
            models = all_models)

plot_figure(ps, variable = "Signal", fig_title =  paste0(path_output, "/Power calc variations/Figure - Signal"),
            models = all_models)

plot_figure(ps, variable = "Magnitude", fig_title =  paste0(path_output, "/Power calc variations/Figure - Magnitude"),
            models = all_models)
