# Helper script for reading the zip files with the results

# Loading required libraries
library(zip)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# Defining the functions to extract simulation results from the zip files
get.figure.data.eval = function(datapath, only_published = T) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  df = ldply(zips, get.data.from.zip.eval, only_published = only_published, .progress = "text")
  df
}

get.figure.data.rep = function(datapath) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  df = ldply(zips, get.data.from.zip.rep, .progress = "text")
  df
}

get.data.from.zip.eval = function(zipfile, only_published = T) {
  
  # There's probably a better way of doing this, but I'll figure it out later,
  # performance is not critical at this point
  source("eval.r", local = T)
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  estimates.df = read.table(
    file = paste(tmpdir, "/estimates.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  
  eval.df = make.evaluation.tests() %>% filter(Published.Only == only_published)
  
  eval.df = dcast(eval.df, . ~ Measure + Statistic, value.var = "Value")
  param.df = dcast(param.df, . ~ Parameter, value.var = "Value")
  
  d = cbind(param.df, eval.df)
  
  return (d)
}

get.data.from.zip.rep = function(zipfile) {
  
  # There's probably a better way of doing this, but I'll figure it out later,
  # performance is not critical at this point
  source("eval.r", local = T)
  source("rep.r", local = T)
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  
  param.df = dcast(param.df, . ~ Parameter, value.var = "Value")
  
  replications.df = data.table(
    read.table(file = paste(tmpdir, "/replications.csv", sep = ""),
               sep = ";", stringsAsFactors = F, header = T))
  
  rep.eval.df = make.rep.evaluation.tests(param.df$min.effect.of.interest)
  
  rep.eval.df = cbind(
    dcast(rep.eval.df %>% filter(!is.na(Type)) %>% select(-LongType),
          RepSet ~ Type + name + N, value.var = "value"),
    dcast(rep.eval.df %>% filter(is.na(Type)) %>% select(-Type, -LongType),
          RepSet ~ name + N, value.var = "value") %>% select(-RepSet)
  )
  
  d = cbind(param.df, rep.eval.df)
  
}