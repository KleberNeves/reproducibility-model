library(shiny)
library(ggplot2)
library(broom)
library(glue)
library(plyr)
library(dplyr)
library(reshape2)
library(metafor)
library(zip)
library(pwr)
library(data.table)

possible.outcomes = c(
  "False Positives",
  "False Negatives",
  "Signal Errors",
  "Exaggeration Factor",
  "True Positives",
  "True Negatives",
  "Potential Discoveries",
  "Discovery Efficiency",
  "Discovery Effectiveness",
  "Reproducibility Rate (BRI)",
  "Reproducibility Rate (BRI) All Papers",
  "Reproducibility Rate (SSS)",
  "Reproducibility Rate (SSS) All Papers"
)

effects.df = data.frame()
estimates.df = data.frame()
eval.df = data.frame()
param.df = data.frame()
master.rep.df = data.frame()

dichotomous.input = list(
  sdA = 0,
  weightB = 0.5,
  meanB = 1,
  sdB = 0,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  alpha.threshold = 0.05,
  n.scientists = 1,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  rep.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = ""
)

two.peaks.input = list(
  sdA = 0.1,
  weightB = 0.5,
  meanB = 1,
  sdB = 0.1,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  alpha.threshold = 0.05,
  n.scientists = 1,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  rep.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = ""
)

single.dist.input = list(
  sdA = 1,
  weightB = 0,
  meanB = 0,
  sdB = 0,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  alpha.threshold = 0.05,
  n.scientists = 1,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  rep.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = ""
)

double.dist.input = list(
  sdA = 0.1,
  weightB = 0.5,
  meanB = 0,
  sdB = 1,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  alpha.threshold = 0.05,
  n.scientists = 1,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  rep.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = ""
)

preset.scenarios = c(
  "Optimistic: low variance and error, large effects, large samples, no publication bias",
  "Optimistic, but with publication bias",
  "Small effects, large sample sizes",
  "Large effects, small sample sizes",
  "Large measurement error",
  "Large interlab variation",
  "Small sample size, small effects",
  "Pessimistic: high variance and error, small sample size, small effects",
  "Small (for quick testing)"
)

shiny_running = TRUE

# Model code
source("model.r")