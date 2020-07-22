# Loading required libraries
library(zip)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)

# Defining the functions to extract simulation results from the zip files
make.diagnostic.plots = function(datapath) {
  zips = paste(datapath, list.files(datapath, ".zip$"), sep ="")
  ps = llply(zips, make.diagnostic.plot, .progress = "text")
  ps
}

make.diagnostic.plot = function(zipfile) {
  
  tmpdir = tempdir()
  unzip(zipfile = zipfile, exdir = tmpdir)
  
  param.df = read.table(
    file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  param.df = dcast(param.df, . ~ Parameter, value.var = "Value")
  
  estimates.df = read.table(
    file = paste(tmpdir, "/estimates.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
  
  D = estimates.df %>%
    select(Effect.Index, Real.Effect.Size, Estimated.Effect.Size,
           p.value, Alpha, Biased) %>%
    mutate(Significant = p.value <= Alpha, Interesting = abs(Real.Effect.Size) > 0.5) %>%
    rename(`Real Effect` = Real.Effect.Size, Estimate = Estimated.Effect.Size) %>%
    select(-Alpha, -p.value) %>%
    melt(id.vars = c("Effect.Index", "Significant", "Interesting", "Biased"))
  
  D$Biased = ifelse(D$Biased, "Biased", "Unbiased")
  D$Interesting = ifelse(D$Interesting, "'True' effects (> 0.5)", "'Null' effects (< 0.5)")
  
  efs = unique(D$Effect.Index)
  efsample = sample(efs, min(1000, length(efs)))
  D = D %>% filter(Effect.Index %in% efsample)
  
  ptitle = paste(
    param.df$scenarioName,
    ", Power = ", param.df$typical.power,
    ", Bias = ", param.df$bias.level,
    ", Prev = ", param.df$`% Above minimum`
  )
  
  p = ggplot(D, aes(x = variable, y = value, group = Effect.Index, color = Significant)) +
      geom_line(alpha = 0.15) + coord_flip() + facet_grid(Interesting~Biased) +
      labs(x = "Effect Size", y = "", title = ptitle)
  
  pw = ifelse(length(unique(D$Biased)) > 1, 11, 7)
  ggsave(paste0(str_replace_all(ptitle, "[.]", "_"), ".png"), p, width = pw, height = 5)
  
  p
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_dir = "/home/kleber/Dropbox/Scientific Research/Projects/Modelo Reprodutibilidade/Results/0720_Paper_Figures/"
diag.plots = make.diagnostic.plots(data_dir)
