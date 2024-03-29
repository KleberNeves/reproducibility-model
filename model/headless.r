# Function that runs a simulation and saves the results
# a.input - a list with the input parameters
# saveFilename - the name of the zip file where results will be saved
hlrun = function(a.input, saveFilename) {
  print(a.input)
  shiny_running <<- FALSE
  run.simulation(a.input)
  dir.create(results_folder)
  
  # Zipping to export results
  fs = c()
  tmpdir = tempdir()
  fn = paste(tmpdir,"/estimates.csv", sep = ""); fs = c(fs, fn)
  write.table(x = estimates.df, file = fn, sep = ";", row.names = F)
  fn = paste(tmpdir,"/eval.csv", sep = ""); fs = c(fs, fn)
  write.table(x = eval.df, file = fn, sep = ";", row.names = F)
  fn = paste(tmpdir,"/pars.csv", sep = ""); fs = c(fs, fn)
  write.table(x = param.df, file = fn, sep = ";", row.names = F)
  
  fn = paste0(tmpdir, "/input.RData"); fs = c(fs, fn)
  save(a.input, file = fn)
  
  if (a.input$calc.repro) {
    fn = paste(tmpdir,"/replications.csv", sep = ""); fs = c(fs, fn)
    write.table(x = replications.df, file = fn, sep = ";", row.names = F)
  }
  
  uldist = data.frame(x = sample.from.dist(a.input$sdA, a.input$weightB,
                                           a.input$meanB, a.input$sdB, k = 100000
  ))
  
  p = ggplot(uldist, aes(x = x)) + geom_density() +
    geom_vline(xintercept = a.input$min.effect.of.interest,
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = -a.input$min.effect.of.interest,
               linetype = "dashed", color = "blue") +
    theme_minimal()
  fn = paste(tmpdir,"/effects_distribution.png", sep = ""); fs = c(fs, fn)
  ggsave(plot = p, file = fn, dpi = 150, width = 7, height = 4)
  
  zipr(zipfile = paste(results_folder, saveFilename, ".zip", sep=""), files = fs)
  
  t = "Saved!"; print(t)
  
  return (0)
}

# Runs a combination of parameters given by a dataframe, substituting the values in the input parameteres list
# Can use a global counter or a local counter (column i in the dataframe) - the counter is used to name the zip files with the results
hlrun_comb = function(comb, global_counter = T) {
  if (global_counter) {
    # Assumes i is a counter existing outside its scope (need to change this)
    i <<- i + 1
    li = i
  } else {
    li = comb$i
    comb$i = NULL
  }
  
  k = baseline.input
  for (n in colnames(comb)) {
    k[n] = as.list(comb)[n]
  }
  k
  
  hlrun(k, paste("Scenario", li))
}

# Saves the information describing that batch of simulations
save.folder.sim.info = function (sim.folder) {
  dir.create(sim.folder)
  
  fn = paste0(sim.folder, "/baseline_input.RData")
  if (!file.exists(fn)) save(baseline.input, file = fn)
  
  fn = paste0(sim.folder, "/run_list.RData")
  if (!file.exists(fn)) save(run.list, file = fn)
  
  fn = paste0(sim.folder, "/model_path.txt")
  if (!file.exists(fn)) {
    f = file(fn)
    writeLines(c(dirname(rstudioapi::getActiveDocumentContext()$path),
                 paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/global.r")), f)
    close(f)
  }
}

# Builds a data frame with all possible combinations of the run.list parameters
# e.g. if given run.list = list(p1 = c(a,b), p2 = c(1,2,3)), it will produce:
# p1 p2
#  a  1
#  a  2
#  a  3
#  b  1
#  b  2
#  b  3
# with all the other parameters equal to the baseline.input
make.full.sweep.df = function (run.list, run_indicator = F) {
  d = expand.grid(run.list, stringsAsFactors = F)
  if (run_indicator) {
    d$i = 1:nrow(d)
    d$has_run = F
  }
  d = post_process_sweep_df(d, run_indicator)
  d
}

# Builds a data frame with combinations of parameters, varying one parameter at a time
# e.g. if given run.list = list(p1 = c(a,b), p2 = c(1,2,3), p3 = c(w,z)), it will produce:
# p1 p2 p3
#  a  1  w
#  b  1  w
#  a  2  w
#  a  3  w
#  a  1  z
#  a  1  z
# with all the other parameters equal to the baseline.input
make.sweep.df = function (run.list, run_indicator = F) {
  d = do.call(rbind, lapply(names(run.list), get_par_sweep))
  d = d[!duplicated(d),]
  if (!is.data.frame(d)) {
    d = as.data.frame(d)
    colnames(d) = names(run.list)
  }
  d = post_process_sweep_df(d, run_indicator)
  d
}

# It has helpers for setting the underlying distributions using scenarioName as a parameter. Specifying any of the distribution parameters ("sdA","sdB","meanB","weightB") in the sweep overrides this helper.
post_process_sweep_df = function (d, run_indicator) {
  if (any(colnames(d) == "scenarioName") & !any(colnames(d) %in% c("sdA","sdB","meanB","weightB"))) {
    if (all(d$scenarioName %in% c("Single Normal", "Two Peaks"))) {
      d$sdA = ifelse(
        d$scenarioName == "Single Normal", 1,
        ifelse(
          d$scenarioName == "Two Peaks", 0.1,
          NA)
      )
      d$sdB = ifelse(
        d$scenarioName == "Single Normal", 0,
        ifelse(
          d$scenarioName == "Two Peaks", 0.1,
          NA)
      )
      d$meanB = ifelse(
        d$scenarioName == "Single Normal", 0,
        ifelse(
          d$scenarioName == "Two Peaks", 1,
          NA)
      )
      d$weightB = ifelse(
        d$scenarioName == "Single Normal", 0,
        ifelse(
          d$scenarioName == "Two Peaks", 0.5,
          NA)
      )
    }
  }
  
  if (run_indicator) {
    d$i = 1:nrow(d)
    d$has_run = F
  }
  
  d
}

# Assumes the existence of a run.df data frame with the combinations
run_sweep_df = function () {
  if (!exists("run.df")) stop("run.df not found.")
  walk(1:nrow(run.df), function (i) {
    if (!run.df[i, "has_run"]) {
      hlr = hlrun_comb(run.df[i, which(colnames(run.df) != "has_run")], global_counter = F)
      if (hlr == 0) run.df[i, "has_run"] <<- T
    }
  })
}

# Generates combinations for a single parameter sweep data frame, all others parameters kept at baseline
get_par_sweep = function(a.par) {
  sweep = run.list
  
  # sdAB is a helper to specify both sdA and sdB to the same values, simultaneously varying
  par_name = a.par
  if (par_name == "sdAB") {
    par_name = c("sdA","sdB")
  }
  if ("sdAB" %in% names(sweep)) {
    sweep["sdA"] = sweep["sdAB"]
    sweep["sdB"] = sweep["sdAB"]
    sweep["sdAB"] = NULL
  }
  
  for (i in 1:length(sweep)) {
    this.par = names(sweep)[i]
    if (!(this.par %in% par_name)) {
      sweep[i] = baseline.input[this.par]
    }
  }
  sweep = compact(sweep)
  r = expand.grid(sweep)
  
  # sdAB helper
  if (a.par == "sdAB") { r = r %>% filter(sdA == sdB) }
  
  r
}
