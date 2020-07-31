library(data.table)
library(metafor)

# Returns a data frame with the results of replications
perform.replications = function(input, rep.power = -1) {
  # Filters the published estimates
  rep.ests = estimates.df[Published == T & p.value <= Alpha]
  
  # Number of replications per experiment
  n.reps = 3
  
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
  browser()
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
    dcast(Effect.Index + Type + LongType + RepSet ~ Measure, value.var = "Value")
  sapply(5:13, function(x) { cast.df[[x]] <<- as.logical(cast.df[[x]]) })
  
  # overall replication rate, PPV/precision, recall
  rep.rates = cast.df[, .(ReproRate = mean(Success),
                          Specificity = sum(TN) / (sum(TN) + sum(FP)),
                          Sensitivity = sum(TP) / (sum(TP) + sum(FN)),
                          SpecificityBias = sum(TNBias) / (sum(TNBias) + sum(FPBias)),
                          SensitivityBias = sum(TPBias) / (sum(TPBias) + sum(FNBias))),
                      by = .(RepSet, Type, LongType)]
  
  # Prevalence of the literature (only of the target)
  cast.df = target.df %>% filter(Measure == "Is.Real") %>% select(-Type, -LongType) %>%
    dcast(Effect.Index + RepSet ~ Measure, value.var = "Value")
  cast.df$Is.Real = as.logical(cast.df$Is.Real)
  other.rates = cast.df[, .(Prev_Sample = mean(Is.Real)), by = .(RepSet)]
  
  cast.df = target.df %>% filter(Measure == "Is.Biased") %>% select(-Type, -LongType) %>%
    dcast(Effect.Index + RepSet ~ Measure, value.var = "Value")
  cast.df$Is.Biased = as.logical(cast.df$Is.Biased)
  other.rates.bias = cast.df[, .(Prev_SampleBias = mean(!Is.Biased)), by = .(RepSet)]
  
  # Prevalence of the literature (whole literature)
  cast.df = rates.df %>% filter(Measure == "Is.Real") %>% select(-Type, -LongType) %>%
    dcast(Effect.Index + RepSet ~ Measure, value.var = "Value")
  cast.df$Is.Real = as.logical(cast.df$Is.Real)
  other.rates.whole = cast.df[, .(Prev_Whole = mean(Is.Real)), by = .(RepSet)]
  
  cast.df = rates.df %>% filter(Measure == "Is.Biased") %>% select(-Type, -LongType) %>%
    dcast(Effect.Index + RepSet ~ Measure, value.var = "Value")
  cast.df$Is.Biased = as.logical(cast.df$Is.Biased)
  other.rates.whole.bias = cast.df[, .(Prev_WholeBias = mean(!Is.Biased)), by = .(RepSet)]
  
  # General (criteria-independent)
  cast.df = target.df %>% filter(is.na(Type) & !(Measure %in% c("Is.Biased", "Is.Real"))) %>%
    select(-Type, -LongType) %>%
    dcast(Effect.Index + RepSet ~ Measure, value.var = "Value")
  sapply(3:4, function(x) { cast.df[[x]] <<- as.numeric(cast.df[[x]]) })
  sapply(5:6, function(x) { cast.df[[x]] <<- as.logical(cast.df[[x]]) })
  
  # exaggeration and signal error rate
  error.rates = cast.df[,
                        .(Exaggeration_MA_x_Original = median(`Exaggeration (MA x Original)`, na.rm = T),
                          Exaggeration_MA_x_Real = median(`Exaggeration (MA x Real)`, na.rm = T),
                          Signal_Error_MA_x_Original = mean(`Signal Error (MA x Original)`, na.rm = T),
                          Signal_Error_MA_x_Real = mean(`Signal Error (MA x Real)`, na.rm = T)),
                        by = .(RepSet)]
  
  # Make long, bind and return
  error.rates = melt(error.rates, id.vars = "RepSet")
  other.rates = melt(other.rates, id.vars = "RepSet")
  other.rates.bias = melt(other.rates.bias, id.vars = "RepSet")
  other.rates.whole = melt(other.rates.whole, id.vars = "RepSet")
  other.rates.whole.bias = melt(other.rates.whole.bias, id.vars = "RepSet")
  rep.rates = melt(rep.rates, id.vars = c("RepSet", "Type", "LongType"))
  
  other.rates = rbind(error.rates, other.rates, other.rates.whole, other.rates.bias, other.rates.whole.bias)
  other.rates$Type = NA
  other.rates$LongType = NA
  
  final.rates = rbind(rep.rates, other.rates)
  final.rates$N = ifelse(n.sample == -1, "All", as.character(n.sample))
  return (final.rates)
}

# For each experiment, calculates reproducibility measures
calc.rep.measures = function(types = c("Orig-in-MA-PI", "MA-SSS", "VOTE-SSS", "Orig-in-MA-CI", "MA-in-Orig-CI"), min.effect.of.interest) {
  # For each set of experiments, calculates each reproducibility measure in types
  rates.df = replications.df[, evaluate.exp.rep(.SD, types = types,
                                                min.effect.of.interest = min.effect.of.interest),
                             by = .(Effect.Index, RepSet)]
  return (rates.df)
}

# Computes many types of reproducibility measures from a set of replications
evaluate.exp.rep = function (rep.exps, types, min.effect.of.interest) {
  
  MA = with(rep.exps, run.ma(MeanControl, SDControl, Sample.Size,
                             MeanTreated, SDTreated, Sample.Size))
  
  result = ldply(types, reproducibility.success, rep.exps = rep.exps, MA = MA)
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
  result = melt(result, id.vars = c("Type", "LongType"))
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
                 
                 # Exaggeration (MA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (MA x Original)",
                            Value = ifelse(MA$m$beta[[1]] / original.estimate <= 0,
                                           NA, MA$m$beta[[1]] / original.estimate)),
                 # Exaggeration (MA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Exaggeration (MA x Real)",
                            Value = ifelse(MA$m$beta[[1]] / real.effect <= 0,
                                           NA, MA$m$beta[[1]] / real.effect)),
                 # Signal (MA x Original)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (MA x Original)",
                            Value = MA$m$beta[[1]] / original.estimate <= 0 &
                              MA$m$pval < 0.05),
                 # Signal (MA x Real)
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Signal Error (MA x Real)",
                            Value = MA$m$beta[[1]] / real.effect <= 0 &
                              MA$m$pval < 0.05)
  )
  
  result
}

# Computes success or failure in a replication according to a give criterium
reproducibility.success = function (type, rep.exps, MA) {
  
  if (type == "Orig-in-MA-PI") {
    longtype = "original estimate is within the PI of the meta analysis"
    value = rep.exps$Original.Effect.Size[1] > MA$pred$cr.lb &
      rep.exps$Original.Effect.Size[1] < MA$pred$cr.ub
  } else if (type == "Orig-in-MA-CI") {
    longtype = "original estimate is within the CI of the meta analysis"
    value = rep.exps$Original.Effect.Size[1] > MA$pred$ci.lb &
      rep.exps$Original.Effect.Size[1] < MA$pred$ci.ub
  } else if (type == "MA-SSS") {
    longtype = "meta analysis is significant and in the same sense as the original"
    value = rep.exps$Original.Effect.Size[1] / MA$m$beta[[1]] > 0 & MA$m$pval < 0.05
  } else if (type == "VOTE-SSS") {
    longtype = "majority of individual replications are significant and in the same sense"
    value = mean(rep.exps$p.value < 0.05) >= 0.5
  } else if (type == "MA-in-Orig-CI") {
    longtype = "meta analysis point estimate is within the CI of the original"
    value = MA$m$beta[[1]] > rep.exps$Original.CI.low[1] & MA$m$beta[[1]] < rep.exps$Original.CI.high[1]
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
run.ma = function(mean_control, sd_control, n_control, mean_treated, sd_treated, n_treated) {
  ess = escalc(measure = "MD", m1i = as.numeric(mean_treated), 
               m2i = as.numeric(mean_control), sd1i = as.numeric(sd_treated), 
               sd2i = as.numeric(sd_control), n1i = as.numeric(n_treated), 
               n2i = as.numeric(n_control))
  tryCatch({
    m = rma(yi = yi, vi = vi, data = ess, measure = "MD", method = "REML",
            control = list(maxiter=1000))
    pred = predict.rma(m, level = 0.95, digits = 3)
    return (list(pred = pred, m = m))
  }, error = function(e) {
    message(e)
    return (NULL)
  })
}