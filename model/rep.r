library(data.table)
library(metafor)
library(tidyr)

# Returns a data frame with the results of replications
perform.replications = function(input, rep.power = -1) {

  # Filters the published estimates
  rep.ests = estimates.df[Published == T & p.value <= Alpha]
  
  # Number of replications per experiment
  n.reps = 3
  if (!is.null(input$repro.exps)) {
    n.reps = input$repro.exps
  }
  
  # Calculate N for the desired power for the original estimate (not the real effect)
  if (rep.power > 0) {
    calc.n = function (eff, sd, wanted.pwr, alpha) {
      pw = tryCatch({
        power.t.test(delta = eff, sd = sd, sig.level = alpha, power = wanted.pwr)$n
      }, error = function(e) {
        # Error means N < 2, i.e. very large effect compared to the SEM
        2
      })
      return (ceiling(pw))
    }
    
    rep.ests$rep.sample.size = mapply(calc.n,
                                      rep.ests$Estimated.Effect.Size,
                                      rep.ests$Estimated.Pooled.SD,
                                      rep.power, input$alpha.threshold)
  } else {
    # Setting rep.power = -1 means: "use the same sample size as the original"
    rep.ests$rep.sample.size = input$typical.sample.size
  }
  
  rep.ests$Original.Effect.Size = rep.ests$Estimated.Effect.Size
  rep.ests$Original.p.value = rep.ests$p.value
  
  # Function to perform a single replication:
  #    receives an effect index, and performs an experiment, adding it to the rep.df
  replicate.exp = function(effect.index, rep.df, separate.reps) {
    rep.input = sanitize_shiny_input(input)
    rep.input$typical.sample.size = rep.ests[Effect.Index == effect.index, rep.sample.size]
    
    # Runs separate repro repeats and saves that information so that we can estimate uncertainty later
    exps = lapply(1:separate.reps, function (rep_i, effect.index, rep.input) {
      exp = perform.experiment(effect.index, rep.input)
      exp$RepSet = rep_i
      exp
    }, effect.index, rep.input)
    
    rbindlist(exps)
  }

  # Runs a number of replications for each effect 
  rep.df = lapply(
    rep(rep.ests$Effect.Index, n.reps),
    replicate.exp, rep.df, separate.reps = input$repro.repeats)
  rep.df = rbindlist(rep.df)
  
  rep.df = merge(rep.df,
                 rep.ests[, .(Effect.Index, rep.sample.size,
                              Original.Effect.Size, CI.low, CI.high, Biased)],
                 by = "Effect.Index")
  
  setnames(rep.df,
           c("CI.low.y","CI.high.y","Biased.y"), c("Original.CI.low", "Original.CI.high","Biased"))
  
  return (rep.df)
}

# Goes over the replication data frame with evaluations and calculates overall summaries for
#   each of the measures present there.
reproducibility.rate = function (rates.df, n.sample = -1) {
  
  # Get a sample if required
  if (n.sample > 0) {
    eff.sample = sample(x = unique(rates.df$Effect.Index), size = n.sample, replace = F)
    target.df = rates.df[Effect.Index %in% eff.sample]
  } else {
    target.df = rates.df
  }
  
  # Summarises replication results (mean, min, max)
  # For each criteria
  cast.df = target.df %>% filter(!is.na(Type)) %>%
    pivot_wider(id_cols = c("Effect.Index", "Type", "LongType", "RepSet"),
                names_from = "Measure", values_from = "Value")
  sapply(5:13, function(x) { cast.df[[x]] <<- as.logical(cast.df[[x]]) })
  
  # overall replication rate, PPV/precision, recall
  cast.df = data.table(cast.df)
  rep.rates = cast.df[, .(ReproRate = mean(Success),
                          Specificity = sum(TN) / (sum(TN) + sum(FP)),
                          Sensitivity = sum(TP) / (sum(TP) + sum(FN)),
                          SpecificityBias = sum(TNBias) / (sum(TNBias) + sum(FPBias)),
                          SensitivityBias = sum(TPBias) / (sum(TPBias) + sum(FNBias))),
                      by = .(RepSet, Type, LongType)]
  
  # Prevalence of the literature (only of the target)
  cast.df = target.df %>% filter(Measure == "Is.Real") %>% select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  cast.df$Is.Real = as.logical(cast.df$Is.Real)
  
  cast.df = data.table(cast.df)
  other.rates = cast.df[, .(Prev_Sample = mean(Is.Real)), by = .(RepSet)]
  
  cast.df = target.df %>% filter(Measure == "Is.Biased") %>% select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  cast.df$Is.Biased = as.logical(cast.df$Is.Biased)
  
  cast.df = data.table(cast.df)
  other.rates.bias = cast.df[, .(Prev_SampleBias = mean(!Is.Biased)), by = .(RepSet)]
  
  # Prevalence of the literature (whole literature)
  cast.df = rates.df %>% filter(Measure == "Is.Real") %>% select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  cast.df$Is.Real = as.logical(cast.df$Is.Real)
  
  cast.df = data.table(cast.df)
  other.rates.whole = cast.df[, .(Prev_Whole = mean(Is.Real)), by = .(RepSet)]
  
  cast.df = rates.df %>% filter(Measure == "Is.Biased") %>% select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  cast.df$Is.Biased = as.logical(cast.df$Is.Biased)
  
  cast.df = data.table(cast.df)
  other.rates.whole.bias = cast.df[, .(Prev_WholeBias = mean(!Is.Biased)), by = .(RepSet)]
  
  # General (criteria-independent)
  cast.df = target.df %>% filter(is.na(Type) & !(Measure %in% c("Is.Biased", "Is.Real"))) %>%
    select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  sapply(c(3,4,7,8), function(x) { cast.df[[x]] <<- as.numeric(cast.df[[x]]) })
  sapply(c(5,6,9,10), function(x) { cast.df[[x]] <<- as.logical(cast.df[[x]]) })
  
  # exaggeration and signal error rate
  cast.df = data.table(cast.df)
  error.rates = cast.df[,
                        .(Exaggeration_RMA_x_Original = median(`Exaggeration (RMA x Original)`,
                                                               na.rm = T),
                          Exaggeration_RMA_x_Real = median(`Exaggeration (RMA x Real)`,
                                                           na.rm = T),
                          Signal_Error_RMA_x_Original = mean(`Signal Error (RMA x Original)`,
                                                             na.rm = T),
                          Signal_Error_RMA_x_Real = mean(`Signal Error (RMA x Real)`,
                                                         na.rm = T),
                          Exaggeration_FMA_x_Original = median(`Exaggeration (FMA x Original)`,
                                                               na.rm = T),
                          Exaggeration_FMA_x_Real = median(`Exaggeration (FMA x Real)`,
                                                           na.rm = T),
                          Signal_Error_FMA_x_Original = mean(`Signal Error (FMA x Original)`,
                                                             na.rm = T),
                          Signal_Error_FMA_x_Real = mean(`Signal Error (FMA x Real)`,
                                                         na.rm = T)),
                        by = .(RepSet)]
  
  # Make long, bind and return
  
  error.rates = pivot_longer(error.rates, cols = -"RepSet")
  other.rates = pivot_longer(other.rates, cols = -"RepSet")
  other.rates.bias = pivot_longer(other.rates.bias, cols = -"RepSet")
  other.rates.whole = pivot_longer(other.rates.whole, cols = -"RepSet")
  other.rates.whole.bias = pivot_longer(other.rates.whole.bias, cols = -"RepSet")
  rep.rates = pivot_longer(rep.rates, cols = -c("RepSet", "Type", "LongType"))
  
  other.rates = rbind(error.rates, other.rates, other.rates.whole, other.rates.bias, other.rates.whole.bias)
  other.rates$Type = NA
  other.rates$LongType = NA
  
  final.rates = rbind(rep.rates, other.rates)
  final.rates$N = ifelse(n.sample == -1, "All", as.character(n.sample))
  return (final.rates)
}

# For each experiment, calculates reproducibility measures
calc.rep.measures = function(types = c("Orig-in-RMA-PI", "RMA-SSS", "Orig-in-RMA-CI", "RMA-in-Orig-CI", "FMA-SSS", "Orig-in-FMA-CI", "FMA-in-Orig-CI", "VOTE-SSS"), min.effect.of.interest) {
  # For each set of experiments, calculates each reproducibility measure in types
  rates.df = replications.df[, evaluate.exp.rep(.SD, types = types,
                                                min.effect.of.interest = min.effect.of.interest),
                             by = .(Effect.Index, RepSet)]
  return (rates.df)
}

# Computes many types of reproducibility measures from a set of replications
evaluate.exp.rep = function (rep.exps, types, min.effect.of.interest) {
  
  RMA = with(rep.exps, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size, type = "RE"))
  
  FMA = with(rep.exps, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size, type = "FE"))
  
  result = ldply(types, reproducibility.success, rep.exps = rep.exps,
                 RMA = RMA, FMA = FMA)
  original.estimate = rep.exps$Original.Effect.Size[1]
  real.effect = rep.exps$Real.Effect.Size[1]
  biased = rep.exps$Biased[1]
  
  # FP (reproduced but is not real)
  result$FP = result$Success & abs(real.effect) < min.effect.of.interest
  # FN (didn't reproduce but is real)
  result$FN = !result$Success & abs(real.effect) >= min.effect.of.interest
  # TP (reproduced and is real)
  result$TP = result$Success & abs(real.effect) >= min.effect.of.interest
  # TN (didn't reproduce and is not real)
  result$TN = !result$Success & abs(real.effect) < min.effect.of.interest
  
  # FP (reproduced but is not real)
  result$FPBias = result$Success & biased
  # FN (didn't reproduce but is real)
  result$FNBias = !result$Success & !biased
  # TP (reproduced and is real)
  result$TPBias = result$Success & !biased
  # TN (didn't reproduce and is not real)
  result$TNBias = !result$Success & biased
  
  # Cast to long format (seems convoluted, there should be a better way)
  result = pivot_longer(result, cols = -c("Type", "LongType"))
  colnames(result) = c("Type", "LongType", "Measure", "Value")
  result$Value = as.character(result$Value)
  
  result = rbind(result,
                 # Real Effect?
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Is.Real",
                            Value = abs(real.effect) >= min.effect.of.interest),
                 
                 # Is Biased?
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Is.Biased",
                            Value = biased),
                 
                 # Exaggeration (RMA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (RMA x Original)",
                            Value = ifelse(RMA$m$beta[[1]] / original.estimate <= 0,
                                           NA, RMA$m$beta[[1]] / original.estimate)),
                 # Exaggeration (RMA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (RMA x Real)",
                            Value = ifelse(RMA$m$beta[[1]] / real.effect <= 0,
                                           NA, RMA$m$beta[[1]] / real.effect)),
                 # Signal (RMA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (RMA x Original)",
                            Value = RMA$m$beta[[1]] / original.estimate <= 0 &
                              RMA$m$pval < 0.05),
                 # Signal (RMA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (RMA x Real)",
                            Value = RMA$m$beta[[1]] / real.effect <= 0 &
                              RMA$m$pval < 0.05),
                 # Exaggeration (FMA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (FMA x Original)",
                            Value = ifelse(FMA$m$beta[[1]] / original.estimate <= 0,
                                           NA, FMA$m$beta[[1]] / original.estimate)),
                 # Exaggeration (FMA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (FMA x Real)",
                            Value = ifelse(FMA$m$beta[[1]] / real.effect <= 0,
                                           NA, FMA$m$beta[[1]] / real.effect)),
                 # Signal (FMA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (FMA x Original)",
                            Value = FMA$m$beta[[1]] / original.estimate <= 0 &
                              FMA$m$pval < 0.05),
                 # Signal (FMA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (FMA x Real)",
                            Value = FMA$m$beta[[1]] / real.effect <= 0 &
                              FMA$m$pval < 0.05)
  )
  
  result
}

# Computes success or failure in a replication according to a give criterium
reproducibility.success = function (type, rep.exps, RMA, FMA) {
  if (type == "Orig-in-RMA-PI") {
    longtype = "original estimate is within the PI of the RE meta analysis"
    value = rep.exps$Original.Effect.Size[1] > RMA$pred$cr.lb &
      rep.exps$Original.Effect.Size[1] < RMA$pred$cr.ub
  } else if (type == "Orig-in-RMA-CI") {
    longtype = "original estimate is within the CI of the RE meta analysis"
    value = rep.exps$Original.Effect.Size[1] > RMA$pred$ci.lb &
      rep.exps$Original.Effect.Size[1] < RMA$pred$ci.ub
  } else if (type == "RMA-SSS") {
    longtype = "RE meta analysis is significant and in the same sense as the original"
    value = rep.exps$Original.Effect.Size[1] / RMA$m$beta[[1]] > 0 & RMA$m$pval < 0.05
  } else if (type == "RMA-in-Orig-CI") {
    longtype = "RE meta analysis point estimate is within the CI of the original"
    value = RMA$m$beta[[1]] > rep.exps$Original.CI.low[1] & RMA$m$beta[[1]] < rep.exps$Original.CI.high[1]
  } else if (type == "Orig-in-FMA-CI") {
    longtype = "original estimate is within the CI of the FE meta analysis"
    value = rep.exps$Original.Effect.Size[1] > FMA$pred$ci.lb &
      rep.exps$Original.Effect.Size[1] < FMA$pred$ci.ub
  } else if (type == "FMA-SSS") {
    longtype = "FE meta analysis is significant and in the same sense as the original"
    value = rep.exps$Original.Effect.Size[1] / FMA$m$beta[[1]] > 0 & FMA$m$pval < 0.05
  } else if (type == "FMA-in-Orig-CI") {
    longtype = "FE meta analysis point estimate is within the CI of the original"
    value = FMA$m$beta[[1]] > rep.exps$Original.CI.low[1] & FMA$m$beta[[1]] < rep.exps$Original.CI.high[1]
  } else if (type == "VOTE-SSS") {
    longtype = "majority of individual replications are significant and in the same sense"
    value = mean(rep.exps$p.value < 0.05) >= 0.5
  } else {
    longtype = "NOT FOUND"
    value = NA
  }
  
  return (data.frame(
    Type = type,
    LongType = longtype,
    Success = value
  ))
  
}

# Runs and returns a meta analysis given the means, SDs and Ns
run.ma = function(mean_control, sd_control, n_control, mean_treated, sd_treated, n_treated, type = "RE") {
  ess = escalc(measure = "MD", m1i = as.numeric(mean_treated), 
               m2i = as.numeric(mean_control), sd1i = as.numeric(sd_treated), 
               sd2i = as.numeric(sd_control), n1i = as.numeric(n_treated), 
               n2i = as.numeric(n_control))
  tryCatch({
    if (type == "RE") {
      m = rma(yi = yi, vi = vi, data = ess, measure = "MD", method = "REML",
              control = list(maxiter=2000, stepadj=0.5))
      pred = predict.rma(m, level = 0.95, digits = 1)
      return (list(pred = pred, m = m))
    } else if (type == "FE") {
      m = rma(yi = yi, vi = vi, data = ess, measure = "MD", method = "FE",
              control = list(maxiter=2000, stepadj=0.5))
      pred = predict.rma(m, level = 0.95, digits = 1)
      return (list(pred = pred, m = m))
    }
  }, error = function(e) {
    message(e)
    print(ess)
    return (NULL)
  })
}
