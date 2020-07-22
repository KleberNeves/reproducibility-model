# Helper script for reading the zip files with the results

# Loading required libraries
library(zip)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# Defining the functions to extract simulation results from the zip files
get.figure.data.rep = function(datapath) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  df = ldply(zips, get.data.from.zip.rep, .progress = "text")
  df
}

get.data.from.zip.rep = function(zipfile) {
  
  # There's probably a better way of doing this, but I'll figure it out later,
  # performance is not critical at this point
  source("../model/eval.r", local = T)
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  estimates.df = read.table(
    file = paste(tmpdir, "/estimates.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  
  eval.df = make.evaluation.tests() %>% filter(Published.Only)
  
  eval.df = dcast(eval.df, . ~ Measure + Statistic, value.var = "Value")
  param.df = dcast(param.df, . ~ Parameter, value.var = "Value")
  
  tryCatch({
    replications.df = read.table(
      file = paste(tmpdir, "/replications.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
    rep.eval.df = make.rep.evaluation.tests(param.df$min.effect.of.interest)
    rep.eval.df = dcast(rep.eval.df, . ~ Measure + Statistic + Type, value.var = "Value")
    d = list(eval = cbind(param.df, eval.df), rep = cbind(param.df, rep.eval.df))
  }, error = function (e) {
    print("No reproducibility data.")
    d = list(eval = cbind(param.df, eval.df))
  })
  return (d)
}