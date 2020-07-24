setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("helper.r")

# Actually getting the data
data_dir = "../../../Results/0720_Paper_Table/"
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

D = figdata
D$UnderlyingDist = str_remove(str_extract(D$scenarioName, ".+;"),";")
D$Scenario = str_remove(str_extract(D$scenarioName, ";.+"),"; ")
D$Scenario = factor(D$Scenario, levels = c("Adequately powered RCT with little bias and 1:1 pre-study odds", "Confirmatory meta-analysis of good quality RCTs", "Meta-analysis of small inconclusive studies", "Underpowered, but well-performed phase I/II RCT", "Underpowered, poorly performed phase I/II RCT", "Adequately powered exploratory epidemiological study", "Underpowered exploratory epidemiological study", "Discovery-oriented exploratory research with massive testing", "Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"))

# PPV
TB1 = D %>% select(power.label, bias.label, `Positive Predictive Value (Correct Signal)_Rate`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Two Normals", "Single Normal")

TB1P = TB1

# Type S
TB1 = D %>% select(power.label, bias.label, `Signal Error_Rate`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Two Normals", "Single Normal")

TB1S = TB1

# Type M
TB1 = D %>% select(power.label, bias.label, `Exaggeration Factor_Median`, Scenario, UnderlyingDist)
TB1 = TB1 %>% melt(id.vars = c("Scenario", "UnderlyingDist","power.label", "bias.label"))
TB1 = TB1 %>% select(-variable) %>%
  dcast(Scenario + power.label + bias.label ~ UnderlyingDist, value.var = "value")

TB1$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB1$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB1 = TB1 %>% select(Scenario, power.label, Odds, Prev, bias.label, Dichotomous, `Two Peaks`, `Two Normals`, `Single Normal`)

colnames(TB1) = c("Scenario description", "Power", "Odds for “True” Effects", "Prevalence of “True” Effects", "Bias", "Dichotomous", "Two Peaks", "Two Normals", "Single Normal")

TB1M = TB1

# TABLE COMPARING ABOVE MIN WITH PREVALENCE

TB2 = D %>% select(`% Above minimum`, Scenario, UnderlyingDist)
TB2 = TB2 %>% melt(id.vars = c("Scenario", "UnderlyingDist"))

TB2 = TB2 %>% select(-variable) %>%
  dcast(Scenario ~ UnderlyingDist, value.var = "value")

TB2$Odds = c("1:1", "2:1", "1:3", "1:5", "1:5", "1:10", "1:10", "1:1000", "1:1000")
TB2$Prev = c("50%", "66,7%", "25%", "16,7%", "16,7%", "9,1%", "9,1%", "0,1%", "0,1%")

TB2 = TB2 %>% select(Scenario, Odds, Prev, Dichotomous, `Two Peaks`, `Two Normals`, `Single Normal`)

colnames(TB2) = c("Scenario description", "Odds for “True” Effects", "Prevalence of “True” Effects", "Dichotomous", "Two Peaks", "Two Normals", "Single Normal")


# Writing the tables to Excel
library(openxlsx)

# Adds column for spacing
TB1P = TB1P %>% mutate(` ` = "") %>% select(1:5,10,6:9)
TB1S = TB1S %>% mutate(` ` = "") %>% select(1:5,10,6:9)
TB1M = TB1M %>% mutate(` ` = "") %>% select(1:5,10,6:9)
TB2 = TB2 %>% mutate(` ` = "") %>% select(1:3,8,4:7)

wb = createWorkbook()

addWorksheet(wb, "PPV")
sht = 1
writeData(wb, sheet = sht, x = TB1P)
setColWidths(wb, sht, cols = 1:10, widths = c(60,10,10,10,10,3,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:10, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:10, gridExpand = TRUE)

addWorksheet(wb, "Type S")
sht = 2
writeData(wb, sheet = sht, x = TB1S)
setColWidths(wb, sht, cols = 1:10, widths = c(60,10,10,10,10,3,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:10, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:10, gridExpand = TRUE)

addWorksheet(wb, "Type M")
sht = 3
writeData(wb, sheet = sht, x = TB1M)
setColWidths(wb, sht, cols = 1:10, widths = c(60,10,10,10,10,3,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:10, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:10, gridExpand = TRUE)

addWorksheet(wb, "Prevalence")
sht = 4
writeData(wb, sheet = sht, x = TB2)
setColWidths(wb, sht, cols = 1:8, widths = c(60,10,10,3,12,12,12,12))
setRowHeights(wb, sht, rows = 1, heights = 40)
setRowHeights(wb, sht, rows = 2:10, heights = 28)
sty = createStyle(valign = "center", halign = "center", textDecoration = "bold", wrapText = T, border = "Bottom", borderStyle = "medium")
addStyle(wb, sheet = sht, sty, rows = 1, cols = 1:8, gridExpand = TRUE)
sty = createStyle(halign = "center", valign = "center", wrapText = T)
addStyle(wb, sheet = sht, sty, rows = 2:10, cols = 1:8, gridExpand = TRUE)

saveWorkbook(wb, "./Table Comparison Ioannidis.xlsx", overwrite = TRUE)
