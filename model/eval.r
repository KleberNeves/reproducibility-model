# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.evaluation.tests = function(repro = F) {
  df = data.frame()
  
  pb = T; df = rbind(df, minimum.effect.count(pb))
  pb = F; df = rbind(df, minimum.effect.count(pb))
  
  pb = T; df = rbind(df, minimum.effect.count.signal(pb))
  pb = F; df = rbind(df, minimum.effect.count.signal(pb))
  
  pb = T; df = rbind(df, true.negatives(pb))
  pb = F; df = rbind(df, true.negatives(pb))
  
  pb = T; df = rbind(df, typeI.error.rate(pb))
  pb = F; df = rbind(df, typeI.error.rate(pb))
  
  pb = T; df = rbind(df, typeII.error.rate(pb))
  pb = F; df = rbind(df, typeII.error.rate(pb))
  
  pb = T; df = rbind(df, typeS.error.rate(pb))
  pb = F; df = rbind(df, typeS.error.rate(pb))
  
  pb = T; df = rbind(df, exaggeration.factor(pb))
  pb = F; df = rbind(df, exaggeration.factor(pb))
  
  pb = T; df = rbind(df, typeS.error.rate.above.min(pb))
  pb = F; df = rbind(df, typeS.error.rate.above.min(pb))
  
  pb = T; df = rbind(df, exaggeration.factor.above.min(pb))
  pb = F; df = rbind(df, exaggeration.factor.above.min(pb))
  
  pb = T; df = rbind(df, discovered.effect.sizes(pb))
  pb = F; df = rbind(df, discovered.effect.sizes(pb))
  
  pb = T; df = rbind(df, pos.pred.value(pb))
  pb = F; df = rbind(df, pos.pred.value(pb))
  
  pb = T; df = rbind(df, pos.pred.value.signal(pb))
  pb = F; df = rbind(df, pos.pred.value.signal(pb))
  
  pb = T; df = rbind(df, neg.pred.value(pb))
  pb = F; df = rbind(df, neg.pred.value(pb))

  if (repro) {
    
    df$variable = character(nrow(df))
    
    t = "Reproducibility ..."
    print(t); if (shiny_running) { showNotification(t, type = "default") }
    
    # Run RR many times and store the results
    make.reps = function (rep_i) {
      
      print(paste(rep_i, "/", input$repro.repeats, sep = ""))
      
      rdf = data.frame(
        Measure = character(0),
        Value = numeric(0),
        Sample.Size = numeric(0),
        Rep = numeric(0)
      )
      
      master.rep.df <<- perform.replications(input, rep.power = 0.95)
      
      # Computes reproducibility rates
      rdf = rbind(rdf, reproducibility.rate(master.rep.df, "BRI", input, n.sample = 20, measure.name = "Reproducibility Rate (BRI)", rep_i))
      rdf = rbind(rdf, reproducibility.rate(master.rep.df, "BRI", input, n.sample = -1, measure.name = "Reproducibility Rate (BRI) All Papers", rep_i))
      rdf = rbind(rdf, reproducibility.rate(master.rep.df, "SSS", input, n.sample = 20, measure.name = "Reproducibility Rate (SSS)", rep_i))
      rdf = rbind(rdf, reproducibility.rate(master.rep.df, "SSS", input, n.sample = -1, measure.name = "Reproducibility Rate (SSS) All Papers", rep_i))
      # rdf = rbind(rdf, reproducibility.rate(master.rep.df, "ST", input, n.sample = 20, measure.name = "Reproducibility Rate (ST)", rep_i))
      # rdf = rbind(rdf, reproducibility.rate(master.rep.df, "ST", input, n.sample = -1, measure.name = "Reproducibility Rate (ST) All Papers", rep_i))
      
      rdf
    }
    
    rr.df = adply(1:input$repro.repeats, 1, make.reps)
    
    # CALCULATE CIs
    rr.df = rr.df %>% group_by(Measure) %>% summarise(MeanValue = mean(Value, na.rm = T),
                                                      SDValue = sd(Value, na.rm = T),
                                                      SEMValue = SDValue / sqrt(input$repro.repeats),
                                                      CI.side = -SEMValue * qnorm(0.05 / 2),
                                                      CI.low = max(0, MeanValue - CI.side, na.rm = T),
                                                      CI.high = min(MeanValue + CI.side, 1, na.rm = T),
                                                      TypeIErrorRateMean = mean(TypeIErrorRate, na.rm = T),
                                                      TypeIIErrorRateMean = mean(TypeIIErrorRate, na.rm = T),
                                                      TypeSErrorRateMean = mean(TypeSErrorRate, na.rm = T),
                                                      MFactorMean = mean(MFactor, na.rm = T),
                                                      PPVMean = mean(PPV, na.rm = T),
                                                      NPVMean = mean(NPV, na.rm = T)
    )
    
    # ADD CIS TO RESULTs
    rr.df = rr.df %>% select(Measure, MeanValue, CI.low, CI.high, TypeIErrorRateMean, TypeIIErrorRateMean, TypeSErrorRateMean, MFactorMean, PPVMean, NPVMean) %>% melt(id.vars = "Measure") %>% mutate(Published.Effects.Only. = T) %>% select(1,3,4,2)
    colnames(rr.df) = colnames(df)
    
    df = rbind(df, rr.df)
    
  } else { # if not calc.repro, just add the variable column
    df$variable = ""
  }
  
  df
}

prop.sep = function (p, n) {
  return ((p * (1 - p) / n) ^ 0.5)
}

# TYPE I/FALSE POSITIVE ERROR RATE
typeI.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "False Positives",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # False Positives
  ests$isError = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  error.count = nrow(ests %>% filter(isError))
  
  value = error.count / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE II/FALSE NEGATIVE ERROR RATE
typeII.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "False Negatives",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # False Negatives
  ests$isError = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  error.count = nrow(ests %>% filter(isError))
  
  value = error.count / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE S/WRONG SIGNAL ERROR RATE
typeS.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "Signal Error",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value < Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Wrong Signal
  ests$isError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  error.count = nrow(ests %>% filter(isError))
  
  value = error.count / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE M/MAGNITUDE/EXAGGERATION FACTOR
exaggeration.factor = function(published.only = T) {
  
  result = data.frame(Measure = "Exaggeration Factor",
                      Statistic = c("Median", "Q25", "Q75", "SEP", "N"),
                      Value = c(NA,NA,NA,NA,0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value < Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Magnitude Error
  ests$Mfactor = ifelse(
    ests$Real.Effect.Size / ests$Estimated.Effect.Size > 0,
    ests$Estimated.Effect.Size / ests$Real.Effect.Size,
    NA)
  
  result$Value = c(median(ests$Mfactor, na.rm = T),
                   quantile(ests$Mfactor, 0.25, na.rm = T),
                   quantile(ests$Mfactor, 0.75, na.rm = T),
                   sd(ests$Mfactor, na.rm = T),
                   nrow(ests))
  return (result)
}

# TYPE S/WRONG SIGNAL ERROR RATE FOR EFFECTS ABOVE MIN OF INTEREST
typeS.error.rate.above.min = function(published.only = T) {
  
  result = data.frame(Measure = "Signal Error (Effects > Min)",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value < Alpha, abs(Estimated.Effect.Size) >= Min.Interesting.Effect)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Wrong Signal
  ests$isError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  error.count = nrow(ests %>% filter(isError))
  
  value = error.count / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE M/MAGNITUDE/EXAGGERATION FACTOR FOR EFFECTS ABOVE MIN OF INTEREST
exaggeration.factor.above.min = function(published.only = T) {
  
  result = data.frame(Measure = "Exaggeration Factor (Effects > Min)",
                      Statistic = c("Median", "Q25", "Q75", "SEP", "N"),
                      Value = c(NA,NA,NA,NA,0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value < Alpha, abs(Real.Effect.Size) >= Min.Interesting.Effect)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Magnitude Error
  ests$Mfactor = ifelse(
    ests$Real.Effect.Size / ests$Estimated.Effect.Size > 0,
    ests$Estimated.Effect.Size / ests$Real.Effect.Size,
    NA)
  
  result$Value = c(median(ests$Mfactor, na.rm = T),
                   quantile(ests$Mfactor, 0.25, na.rm = T),
                   quantile(ests$Mfactor, 0.75, na.rm = T),
                   sd(ests$Mfactor, na.rm = T),
                   nrow(ests))
  return (result)
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE
pos.pred.value = function(published.only = T) {
  
  result = data.frame(Measure = "Positive Predictive Value",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value <= Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Positive Predictive Value = True Positives / (True Positives + False Positives)
  ests$isTP = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  tps = nrow(ests %>% filter(isTP))
  
  value = tps / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE (CONDITIONAL ON SIGNAL BEING RIGHT)
pos.pred.value.signal = function(published.only = T) {
  
  result = data.frame(Measure = "Positive Predictive Value (Correct Signal)",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value <= Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Wrong Signal
  ests$signalError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  ests = ests %>% filter(!signalError)
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Positive Predictive Value = True Positives / (True Positives + False Positives)
  ests$isTP = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) & (ests$p.value <= ests$Alpha)
  
  tps = nrow(ests %>% filter(isTP))
  
  value = tps / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# NEGATIVE PREDICTIVE VALUE
neg.pred.value = function(published.only = T) {
  
  result = data.frame(Measure = "Negative Predictive Value",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value > Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Negative Predictive Value = True Negatives / (True Negatives + False Negatives)
  ests$isTN = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  tns = nrow(ests %>% filter(isTN))
  
  value = tns / nrow(ests)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TRUE POSITIVES
minimum.effect.count = function(published.only = T) {
  
  result = data.frame(Measure = "True Positives",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # True positives
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery)) / nrow(ests)
  
  value = count
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# TRUE POSITIVES (CONDITIONAL ON SIGNAL BEING RIGHT)
minimum.effect.count.signal = function(published.only = T) {
  
  result = data.frame(Measure = "True Positives (Correct Signal)",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Wrong Signal
  ests$signalError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  ests = ests %>% filter(!signalError)
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # True positives
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery)) / nrow(ests)
  
  value = count
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}


# TRUE NEGATIVES
true.negatives = function(published.only = T) {
  
  result = data.frame(Measure = "True Negatives",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # True negatives
  ests$isNonDiscovery = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  count = nrow(ests %>% filter(isNonDiscovery)) / nrow(ests)
  
  value = count
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# N TRUE EFFECTS
n.real.effects.above.minimum = function(published.only = T) {
  
  result = data.frame(Measure = "Discoveries",
                      Statistic = c("Number", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Real effects above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect)
  
  count = nrow(ests %>% filter(isDiscovery))
  
  value = round(count,0)
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# Number of true positives / total sample size
efficiency = function(published.only = T) {
  
  result = data.frame(Measure = "Efficiency",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # True positives above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery))
  resources = sum(estimates.df$Sample.Size)
  
  value = count / resources
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# Number of true positives / number of potential discoveries
effectiveness = function(published.only = T) {
  
  result = data.frame(Measure = "Effectiveness",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # True positives above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery))
  tried = n.real.effects.above.minimum(published.only)
  
  value = count / tried
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}

# DISCOVERED EFFECTS' SIZES
discovered.effect.sizes = function(published.only = T) {
  
  result = data.frame(Measure = "DiscoveredES",
                      Statistic = c("Median", "Q25", "Q75", "SEP", "N"),
                      Value = c(NA,NA,NA,NA,0),
                      Published.Only = published.only)
  
  ests = estimates.df %>% filter(p.value < Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Size of effects
  ests$Mfactor = abs(ests$Real.Effect.Size)
  
  result$Value = c(median(ests$Mfactor, na.rm = T),
                   quantile(ests$Mfactor, 0.25, na.rm = T),
                   quantile(ests$Mfactor, 0.75, na.rm = T),
                   sd(ests$Mfactor, na.rm = T),
                   nrow(ests))
  
  return (result)
}
