library(shiny)
library(tidyverse)
library(broom)
library(glue)
library(metafor)
library(zip)
library(pwr)
library(data.table)

global_rep_types = c("VOTE_SSS_005", "VOTE_SSS_0005", "FMA_SSS_005", "FMA_SSS_0005", "RMA_SSS_005", "RMA_SSS_0005", "ORIG_IN_RMA_PI", "ORIG_IN_FMA_CI", "REP_IN_ORIG_CI", "CMA_SSS_005", "CMA_SSS_0005", "SMALL_TELESCOPE", "BF_3", "BF_10")

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

# effects.df = data.table()
# estimates.df = data.table()
# eval.df = data.frame()
# param.df = data.frame()
# replications.df = data.table()

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
  repro.sample = 20,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F,
  fixed.prev.mode = "none",
  fixed.prev = 0
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
  repro.sample = 20,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F,
  fixed.prev.mode = "none",
  fixed.prev = 0
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
  repro.sample = 20,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F,
  fixed.prev.mode = "none",
  fixed.prev = 0
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
  repro.sample = 20,
  repro.power = 0.95,
  repro.detect = 0,
  how.sim.ends = "At a given number of published effects",
  sim.end.value = 5000,
  scenarioName = "",
  publish.only.if.large = F,
  fixed.prev.mode = "none",
  fixed.prev = 0
)

shiny_running = TRUE

sanitize_shiny_input = function (a.input) {
  if ("reactivevalues" %in% class(a.input)) {
    inputlist = reactiveValuesToList(a.input)
  } else {
    inputlist = a.input
  }
  standard_input = single.dist.input
  
  missing = names(single.dist.input)[!(names(single.dist.input) %in% names(inputlist))]
  for (x in missing) {
    inputlist[[x]] = standard_input[[x]]
  }
  
  return (inputlist)
}

preset.scenarios = c("",
                     "Two peaks, high power, low bias",
                     "Two peaks, low power, high bias",
                     "Continuous, high power, low bias",
                     "Continuous, low power, high bias"
)

tryCatch({
  # Help text
  overview_help_text = read_file("./help/overview.html")
  parameters_help_text = read_file("./help/parameters.html")
  outcomes_help_text = read_file("./help/outcomes.html")

  # Source model code
  source("model.r")
  source("eval.r")
  source("rep.r")
}, error = function (e) {
  print("Some files not found. If you're sourcing this for analysis, don't worry.")
  # print(e)
})