# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.evaluation.tests = function() {
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
  
  pb = T; df = rbind(df, pos.pred.value(pb))
  pb = F; df = rbind(df, pos.pred.value(pb))
  
  pb = T; df = rbind(df, pos.pred.value.signal(pb))
  pb = F; df = rbind(df, pos.pred.value.signal(pb))
  
  pb = T; df = rbind(df, neg.pred.value(pb))
  pb = F; df = rbind(df, neg.pred.value(pb))

  df
}

prop.sd = function (p, n) {
  return ((n * p * (1 - p) ^ 0.5))
}

# TYPE I/FALSE POSITIVE ERROR RATE
typeI.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "False Positives",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE II/FALSE NEGATIVE ERROR RATE
typeII.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "False Negatives",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE S/WRONG SIGNAL ERROR RATE
typeS.error.rate = function(published.only = T) {
  
  result = data.frame(Measure = "Signal Error",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE M/MAGNITUDE/EXAGGERATION FACTOR
exaggeration.factor = function(published.only = T) {
  
  result = data.frame(Measure = "Exaggeration Factor",
                      Statistic = c("Median", "IQR", "SD", "N"),
                      Value = c(NA,NA,NA,0))
  
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
                   IQR(ests$Mfactor, na.rm = T),
                   sd(ests$Mfactor, na.rm = T),
                   nrow(ests))
  return (result)
}

# TYPE S/WRONG SIGNAL ERROR RATE FOR EFFECTS ABOVE MIN OF INTEREST
typeS.error.rate.above.min = function(published.only = T) {
  
  result = data.frame(Measure = "Signal Error (Effects > Min)",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TYPE M/MAGNITUDE/EXAGGERATION FACTOR FOR EFFECTS ABOVE MIN OF INTEREST
exaggeration.factor.above.min = function(published.only = T) {
  
  result = data.frame(Measure = "Exaggeration Factor (Effects > Min)",
                      Statistic = c("Median", "IQR", "SD", "N"),
                      Value = c(NA,NA,NA,0))
  
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
                   IQR(ests$Mfactor, na.rm = T),
                   sd(ests$Mfactor, na.rm = T),
                   nrow(ests))
  return (result)
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE
pos.pred.value = function(published.only = T) {
  
  result = data.frame(Measure = "Positive Predictive Value",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE (CONDITIONAL ON SIGNAL BEING RIGHT)
pos.pred.value.signal = function(published.only = T) {
  
  result = data.frame(Measure = "Positive Predictive Value (Correct Signal)",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# NEGATIVE PREDICTIVE VALUE
neg.pred.value = function(published.only = T) {
  
  result = data.frame(Measure = "Negative Predictive Value",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TRUE POSITIVES
minimum.effect.count = function(published.only = T) {
  
  result = data.frame(Measure = "True Positives",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# TRUE POSITIVES (CONDITIONAL ON SIGNAL BEING RIGHT)
minimum.effect.count.signal = function(published.only = T) {
  
  result = data.frame(Measure = "True Positives (Correct Signal)",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}


# TRUE NEGATIVES
true.negatives = function(published.only = T) {
  
  result = data.frame(Measure = "True Negatives",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# N TRUE EFFECTS
n.real.effects.above.minimum = function(published.only = T) {
  
  result = data.frame(Measure = "Discoveries",
                      Statistic = c("Number", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# Number of true positives / total sample size
efficiency = function(published.only = T) {
  
  result = data.frame(Measure = "Efficiency",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}

# Number of true positives / number of potential discoveries
effectiveness = function(published.only = T) {
  
  result = data.frame(Measure = "Effectiveness",
                      Statistic = c("Rate", "SD", "N"),
                      Value = c(NA, NA, 0))
  
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
  result$Value = c(value, prop.sd(value, nrow(ests)), nrow(ests))
  return (result)
}