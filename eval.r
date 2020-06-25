# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.evaluation.tests = function(input) {
  df = data.frame()
  
  t = "Error rates and discoveries ..."
  print(t); if (shiny_running) { showNotification(t, type = "default") }
  
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
  pb = F; df = rbind(df, pos.pred.value(pb),)
  
  pb = T; df = rbind(df, pos.pred.value.signal(pb))
  pb = F; df = rbind(df, pos.pred.value.signal(pb))
  
  pb = T; df = rbind(df, neg.pred.value(pb))
  pb = F; df = rbind(df, neg.pred.value(pb))
  
  if (input$calc.repro) {
    
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

# TYPE I/FALSE POSITIVE ERROR RATE
typeI.error.rate = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NA)
  }
  
  # False Positives
  ests$isError = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  error.count = nrow(ests %>% filter(isError))
  return (error.count / nrow(ests))
}

# TYPE II/FALSE NEGATIVE ERROR RATE
typeII.error.rate = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # False Negatives
  ests$isError = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  error.count = nrow(ests %>% filter(isError))
  return (error.count / nrow(ests))
}

# TYPE S/WRONG SIGNAL ERROR RATE
typeS.error.rate = function(published.only = T) {
  ests = estimates.df %>% filter(p.value < Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Wrong Signal
  ests$isError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  error.count = nrow(ests %>% filter(isError))
  return (error.count / nrow(ests))
}


# TYPE M/MAGNITUDE/EXAGGERATION FACTOR
exaggeration.factor = function(published.only = T) {
  ests = estimates.df %>% filter(p.value < Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Magnitude Error
  ests$Mfactor = ifelse(
    ests$Real.Effect.Size / ests$Estimated.Effect.Size > 0,
    ests$Estimated.Effect.Size / ests$Real.Effect.Size,
    NA)
  
  Mfactor = median(ests$Mfactor, na.rm = T)
  return (Mfactor)
}

# TYPE S/WRONG SIGNAL ERROR RATE FOR EFFECTS ABOVE MIN OF INTEREST
typeS.error.rate.above.min = function(published.only = T) {
  ests = estimates.df %>% filter(p.value < Alpha, abs(Estimated.Effect.Size) >= Min.Interesting.Effect)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Wrong Signal
  ests$isError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  error.count = nrow(ests %>% filter(isError))
  return (error.count / nrow(ests))
}

# TYPE M/MAGNITUDE/EXAGGERATION FACTOR FOR EFFECTS ABOVE MIN OF INTEREST
exaggeration.factor.above.min = function(published.only = T) {
  ests = estimates.df %>% filter(p.value < Alpha, abs(Real.Effect.Size) >= Min.Interesting.Effect)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Magnitude Error
  ests$Mfactor = ifelse(
    ests$Real.Effect.Size / ests$Estimated.Effect.Size > 0,
    ests$Estimated.Effect.Size / ests$Real.Effect.Size,
    NA)
  
  Mfactor = median(ests$Mfactor, na.rm = T)
  return (Mfactor)
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE
pos.pred.value = function(published.only = T) {
  ests = estimates.df %>% filter(p.value <= Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Positive Predictive Value = True Positives / (True Positives + False Positives)
  ests$isTP = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  tps = nrow(ests %>% filter(isTP))
  return (tps / nrow(ests))
}

# POSITIVE PREDICTIVE VALUE/1 - FALSE FINDING RATE (CONDITIONAL ON SIGNAL BEING RIGHT)
pos.pred.value.signal = function(published.only = T) {
  ests = estimates.df %>% filter(p.value <= Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Wrong Signal
  ests$signalError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  ests = ests %>% filter(!signalError)
  
  if (nrow(ests) == 0) {
    return (NA)
  }
  
  # Positive Predictive Value = True Positives / (True Positives + False Positives)
  ests$isTP = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) & (ests$p.value <= ests$Alpha)
  
  tps = nrow(ests %>% filter(isTP))
  return (tps / nrow(ests))
}

# NEGATIVE PREDICTIVE VALUE
neg.pred.value = function(published.only = T) {
  ests = estimates.df %>% filter(p.value > Alpha)
  
  if (published.only) {
    ests = ests %>% filter(Published == T)
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Negative Predictive Value = True Negatives / (True Negatives + False Negatives)
  ests$isTN = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  tns = nrow(ests %>% filter(isTN))
  return (tns / nrow(ests))
}

# TRUE POSITIVES
minimum.effect.count = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # True positives
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery)) / nrow(ests)
  return (count)
}

# TRUE POSITIVES (CONDITIONAL ON SIGNAL BEING RIGHT)
minimum.effect.count.signal = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Wrong Signal
  ests$signalError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  ests = ests %>% filter(!signalError)
  
  if (nrow(ests) == 0) {
    return (NA)
  }
  
  # True positives
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery)) / nrow(ests)
  return (count)
}


# TRUE NEGATIVES
true.negatives = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # True negatives
  ests$isNonDiscovery = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) &
    (ests$p.value > ests$Alpha)
  
  count = nrow(ests %>% filter(isNonDiscovery)) / nrow(ests)
  return (count)
}

# N TRUE EFFECTS
n.real.effects.above.minimum = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # Real effects above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect)
  
  count = nrow(ests %>% filter(isDiscovery))
  return (round(count,0))
}

# Number of true positives / total sample size
efficiency = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # True positives above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery))
  resources = sum(estimates.df$Sample.Size)
  
  return (count / resources)
}

# Number of true positives / number of potential discoveries
effectiveness = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NaN)
  }
  
  # True positives above minimum of interest
  ests$isDiscovery = (abs(ests$Real.Effect.Size) >= ests$Min.Interesting.Effect) &
    (ests$p.value <= ests$Alpha)
  
  count = nrow(ests %>% filter(isDiscovery))
  tried = n.real.effects.above.minimum(published.only)
  
  return (count / tried)
}