# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.evaluation.tests = function() {
  df = data.frame()

  pb = T; df = rbind(df, true.positive.rate(pb))
  pb = F; df = rbind(df, true.positive.rate(pb))
  
  pb = T; df = rbind(df, true.positive.rate.signal(pb))
  pb = F; df = rbind(df, true.positive.rate.signal(pb))
  
  pb = T; df = rbind(df, true.negative.rate(pb))
  pb = F; df = rbind(df, true.negative.rate(pb))
  
  pb = T; df = rbind(df, false.positive.rate(pb))
  pb = F; df = rbind(df, false.positive.rate(pb))
  
  pb = T; df = rbind(df, false.negative.rate(pb))
  pb = F; df = rbind(df, false.negative.rate(pb))
  
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
  
  df
}

# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.rep.evaluation.tests = function(min.effect) {
  
  rates.df = calc.rep.measures(min.effect.of.interest = min.effect)
  
  df = data.frame()
  
  df = rbind(df, reproducibility.rate(rates.df, n.sample = -1))
  df = rbind(df, reproducibility.rate(rates.df, n.sample = 20))

  df
}

prop.sep = function (p, n) {
  return ((p * (1 - p) / n) ^ 0.5)
}

# TYPE I/FALSE POSITIVE ERROR RATE
false.positive.rate = function(published.only = T) {
  
  result = data.frame(Measure = "False Positives",
                      Statistic = c("Rate", "SEP", "N"),
                      Value = c(NA, NA, 0),
                      Published.Only = published.only)
  # browser()
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
false.negative.rate = function(published.only = T) {
  
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
  
  if (nrow(ests) == 0) {
    return (result)
  }
  
  # Positive Predictive Value = True Positives / (True Positives + False Positives)
  ests$isTP = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) & (ests$p.value <= ests$Alpha) & !ests$signalError
  
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
true.positive.rate = function(published.only = T) {
  
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
true.positive.rate.signal = function(published.only = T) {
  
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
  
  # True positives
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha) & !ests$signalError
  
  count = nrow(ests %>% filter(isDiscovery)) / nrow(ests)
  
  value = count
  result$Value = c(value, prop.sep(value, nrow(ests)), nrow(ests))
  return (result)
}


# TRUE NEGATIVES
true.negative.rate = function(published.only = T) {
  
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
