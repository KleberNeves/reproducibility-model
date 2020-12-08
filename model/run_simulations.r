  # Script to run simulations without the Shiny interface
  
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  source("model/global.r")
  source("model/headless.r")
  
  ############## RUNS BELOW ################
  n.sims.to.run = 5000
  
  ######### SIMS FOR FIGURES #########
  results_folder = "./Figures"
  
  i = 0
  
  # Set list of parameters to vary in the sweep (the ones that are not included here will remain fixed across sims)
  run.list = list(
    weightB = c(0.1, 0.3, 0.5, 0.7, 0.9),
    bias.level = c(0, 0.2, 0.6),
    typical.power = c(0.2, 0.5, 0.8)
  )
  
  # Dichotomous
  baseline.input = dichotomous.input
  baseline.input$scenarioName = "Dichotomous"
  baseline.input$typical.sample.size = NULL
  baseline.input$sim.end.value = n.sims.to.run
  
  run.df = make.full.sweep.df(run.list)
  d_ply(run.df, colnames(run.df), hlrun_comb)
  
  # Separate Peaks
  baseline.input = dichotomous.input
  baseline.input$scenarioName = "Peaks SD 0.1"
  baseline.input$typical.sample.size = NULL
  baseline.input$sdA = 0.1
  baseline.input$sdB = 0.1
  baseline.input$sim.end.value = n.sims.to.run
  
  run.df = make.full.sweep.df(run.list)
  d_ply(run.df, colnames(run.df), hlrun_comb)
  
  # Overlapping Peaks
  baseline.input = dichotomous.input
  baseline.input$scenarioName = "Peaks SD 0.3"
  baseline.input$typical.sample.size = NULL
  baseline.input$sdA = 0.3
  baseline.input$sdB = 0.3
  baseline.input$sim.end.value = n.sims.to.run
  
  run.df = make.full.sweep.df(run.list)
  d_ply(run.df, colnames(run.df), hlrun_comb)
    
  # Two Normals
  baseline.input = double.dist.input
  baseline.input$scenarioName = "Two Normals"
  baseline.input$typical.sample.size = NULL
  baseline.input$sdA = 0.1
  baseline.input$sdB = 1
  baseline.input$sim.end.value = n.sims.to.run
  
  run.df = make.full.sweep.df(run.list)
  d_ply(run.df, colnames(run.df), hlrun_comb)
  
  # Single Normal
  run.list = list(
    # Sets SD instead of prevalence/weight. Choice of SD is to approximately match the prevalence for the other distributions.
    sdA = c(0.3, 0.5, 0.8, 1.5, 4),
    bias.level = c(0, 0.2, 0.6),
    typical.power = c(0.2, 0.5, 0.8)
  )
  
  baseline.input = single.dist.input
  baseline.input$scenarioName = "Single Normal"
  baseline.input$typical.sample.size = NULL
  baseline.input$sim.end.value = n.sims.to.run
  
  run.df = make.full.sweep.df(run.list)
  d_ply(run.df, colnames(run.df), hlrun_comb)
  
  ######### SIMS FOR SCENARIO TABLE #########
  results_folder = "./Table"
  
  n.sims.to.run = 5000
  
  i = 0
  
  ### Dichotomous
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Adequately powered RCT with little bias and 1:1 pre-study odds"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.1
  k$weightB = 0.5
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Confirmatory meta-analysis of good quality RCTs"
  k$typical.sample.size = NULL
  k$typical.power = 0.95
  k$bias.level = 0.3
  k$weightB = 0.667
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Meta-analysis of small inconclusive studies"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.4
  k$weightB = 0.25
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Underpowered, but well-performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.167
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Underpowered, poorly performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.167
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Adequately powered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.3
  k$weightB = 0.091
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Underpowered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.3
  k$weightB = 0.091
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Discovery-oriented exploratory research with massive testing"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.001
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))

  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Dichotomous; Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.001
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  
  ### Overlapping Peaks
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Adequately powered RCT with little bias and 1:1 pre-study odds"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.1
  k$weightB = 0.475
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Confirmatory meta-analysis of good quality RCTs"
  k$typical.sample.size = NULL
  k$typical.power = 0.95
  k$bias.level = 0.3
  k$weightB = 0.667
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Meta-analysis of small inconclusive studies"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.4
  k$weightB = 0.185
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Underpowered, but well-performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.085
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Underpowered, poorly performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.085
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Adequately powered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.3
  k$weightB = 0.0001
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Underpowered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.3
  k$weightB = 0.0001
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Discovery-oriented exploratory research with massive testing"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.0001
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Overlapping Peaks; Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.0001
  k$sdA = 0.3
  k$sdB = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  
  ### Two Peaks
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Adequately powered RCT with little bias and 1:1 pre-study odds"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.1
  k$weightB = 0.5
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Confirmatory meta-analysis of good quality RCTs"
  k$typical.sample.size = NULL
  k$typical.power = 0.95
  k$bias.level = 0.3
  k$weightB = 0.667
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Meta-analysis of small inconclusive studies"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.4
  k$weightB = 0.25
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Underpowered, but well-performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.167
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Underpowered, poorly performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.167
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Adequately powered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.3
  k$weightB = 0.091
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Underpowered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.3
  k$weightB = 0.091
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Discovery-oriented exploratory research with massive testing"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.001
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = dichotomous.input
  k$scenarioName = "Two Peaks; Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.001
  k$sdA = 0.1
  k$sdB = 0.1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  ### Two Normals
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Adequately powered RCT with little bias and 1:1 pre-study odds"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.1
  k$weightB = 0.81
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Confirmatory meta-analysis of good quality RCTs"
  k$typical.sample.size = NULL
  k$typical.power = 0.95
  k$bias.level = 0.3
  k$weightB = 1
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Meta-analysis of small inconclusive studies"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.4
  k$weightB = 0.41
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Underpowered, but well-performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.27
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Underpowered, poorly performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.27
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Adequately powered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.3
  k$weightB = 0.15
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Underpowered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.3
  k$weightB = 0.15
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Discovery-oriented exploratory research with massive testing"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$weightB = 0.002
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = double.dist.input
  k$scenarioName = "Two Normals; Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$weightB = 0.002
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  
  ### Single Normal
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Adequately powered RCT with little bias and 1:1 pre-study odds"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.1
  k$sdA = 0.73
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Confirmatory meta-analysis of good quality RCTs"
  k$typical.sample.size = NULL
  k$typical.power = 0.95
  k$bias.level = 0.3
  k$sdA = 1.16
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Meta-analysis of small inconclusive studies"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.4
  k$sdA = 0.44
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Underpowered, but well-performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$sdA = 0.36
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Underpowered, poorly performed phase I/II RCT"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$sdA = 0.36
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Adequately powered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.8
  k$bias.level = 0.3
  k$sdA = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Underpowered exploratory epidemiological study"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.3
  k$sdA = 0.3
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Discovery-oriented exploratory research with massive testing"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.8
  k$sdA = 0.15
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
  
  i <<- i + 1
  k = single.dist.input
  k$scenarioName = "Single Normal; Discovery-oriented exploratory research with massive testing, but with more limited bias (more standardized)"
  k$typical.sample.size = NULL
  k$typical.power = 0.2
  k$bias.level = 0.2
  k$sdA = 0.15
  k$sim.end.value = n.sims.to.run
  hlrun(k,paste("Scenario", i))
