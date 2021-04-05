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
}

hlrun_comb = function(comb) {
  # Assumes i is a counter existing outside its scope (need to change this)
  i <<- i + 1
  k = baseline.input
  for (n in colnames(comb)) {
    k[n] = as.list(comb)[n]
  }
  k
  
  hlrun(k, paste("Scenario", i))
}

save.folder.sim.info = function (sim.folder) {
  fn = paste0(sim.folder, "/baseline_input.RData")
  save(baseline.input, file = fn)
  
  fn = paste0(sim.folder, "/run_list.RData")
  save(run.list, file = fn)
  
  fn = paste0(sim.folder, "/model_path.txt")
  f = file(fn)
  writeLines(c(dirname(rstudioapi::getActiveDocumentContext()$path),
               paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/global.r")), f)
  close(f)
}

make.full.sweep.df = function (run.list) {
  d = expand.grid(run.list, stringsAsFactors = F)
  d
}

make.sweep.df = function (run.list) {
  d = do.call(rbind, lapply(names(run.list), get_par_sweep))
  d = d[!duplicated(d),]
  if (!is.data.frame(d)) {
    d = as.data.frame(d)
    colnames(d) = names(run.list)
  }
  d
}

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
