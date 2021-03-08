library(shiny)
library(ggplot2)
library(broom)
library(glue)
library(purrr)
library(dplyr)
library(tidyr)
library(metafor)
library(zip)
library(pwr)
library(data.table)
library(readr)

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

effects.df = data.table()
estimates.df = data.table()
eval.df = data.frame()
param.df = data.frame()
replications.df = data.table()

dichotomous.input = list(
  sdA = 0,
  weightB = 0.5,
  meanB = 1,
  sdB = 0,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  power.calculation = "mean of the true effects",
  alpha.threshold = 0.05,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  repro.exps = 3,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F
)

two.peaks.input = list(
  sdA = 0.1,
  weightB = 0.5,
  meanB = 1,
  sdB = 0.1,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  power.calculation = "mean of the true effects",
  alpha.threshold = 0.05,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  repro.exps = 3,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F
)

single.dist.input = list(
  sdA = 1,
  weightB = 0,
  meanB = 0,
  sdB = 0,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  power.calculation = "mean of the true effects",
  alpha.threshold = 0.05,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  repro.exps = 3,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F
)

double.dist.input = list(
  sdA = 0.1,
  weightB = 0.5,
  meanB = 0,
  sdB = 1,
  min.effect.of.interest = 0.5,
  typical.sample.size = 17,
  typical.power = NULL,
  power.calculation = "mean of the true effects",
  alpha.threshold = 0.05,
  measure.error = 0,
  interlab.var = 0,
  neg.incentive = 0,
  bias.level = 0,
  calc.repro = F,
  repro.repeats = 5,
  repro.exps = 3,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F
)

shiny_running = TRUE

sanitize_shiny_input = function (a.input) {
  if ("reactivevalues" %in% class(a.input)) {
    return (reactiveValuesToList(a.input))
  } else {
    return (a.input)
  }
}

preset.scenarios = c("",
                     "Two peaks, high power, low bias",
                     "Two peaks, low power, high bias",
                     "Continuous, high power, low bias",
                     "Continuous, low power, high bias"
)

# Help text
overview_help_text = read_file("./help/overview.html")
parameters_help_text = read_file("./help/parameters.html")
outcomes_help_text = read_file("./help/outcomes.html")

# Source model code
source("model.r")
source("eval.r")
source("rep.r")