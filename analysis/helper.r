# Helper script for reading the zip files with the results

# Loading required libraries
library(zip)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)

# Defining the functions to extract simulation results from the zip files
get.figure.data.eval = function(datapath, only_published = T) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  total = length(zips); i = 0
  df = map_dfr(zips, function (z) {
    i <<- i + 1
    print(paste0(i, "/", total))
    get.data.from.zip.eval(z, only_published = only_published)
  })
  df
}

get.figure.data.rep = function(datapath) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  total = length(zips); i = 0
  df = map_dfr(zips, function (z) {
    i <<- i + 1
    print(paste0(i, "/", total))
    get.data.from.zip.rep(z)
  })
  df
}

get.data.from.zip.eval = function(zipfile, only_published = T) {
  
  source("./eval.r", local = T)
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste0(tmpdir, "/pars.csv"), sep = ";", stringsAsFactors = F, header = T)
  estimates.df = read.table(
    file = paste0(tmpdir, "/estimates.csv"), sep = ";", stringsAsFactors = F, header = T)
  
  eval.df = make.evaluation.tests() %>% filter(Published.Only == only_published)
  
  eval.df = eval.df %>% pivot_wider(names_from = c(Measure, Statistic), values_from = Value)
  param.df = param.df %>% pivot_wider(names_from = Parameter, values_from = Value)
  
  d = cbind(param.df, eval.df)
  
  return (d)
}

get.data.from.zip.rep = function(zipfile) {
  
  # There's probably a better way of doing this, but I'll figure it out later,
  # performance is not critical at this point
  source("./eval.r", local = T)
  source("./rep.r", local = T)
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  
  param.df = param.df %>% pivot_wider(names_from = Parameter, values_from = Value)
  
  load(paste0(tmpdir, "/input.RData"))
  
  replications.df = data.table(
    read.table(file = paste0(tmpdir, "/replications.csv"),
               sep = ";", stringsAsFactors = F, header = T))
  
  rep.eval.df = make.rep.evaluation.tests(a.input)
  
  rep.eval.df = cbind(
    rep.eval.df %>%
      filter(!is.na(Type)) %>%
      pivot_wider(id_cols = RepSet, names_from = c(Type, name, N, Nprop), values_from = value),
    rep.eval.df %>%
      filter(is.na(Type)) %>%
      select(-Type) %>%
      pivot_wider(id_cols = RepSet, names_from = c(name, N, Nprop), values_from = value)
  )
  
  d = cbind(param.df, rep.eval.df)
  
  d
}