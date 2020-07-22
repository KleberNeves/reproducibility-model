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
    exps = lapply(1:separate.reps, function (rep_i, ...) {
      exp = perform.experiment(...)
      exp$RepSet = rep_i
    }, effect.index, rep.input, -1)
    
    rbindlist(exps)
  }
  
  # Runs a number of replications for each effect 
  rep.df = lapply(
    rep(rep.ests$Effect.Index, n.reps),
    replicate.exp, rep.df, separate.reps = input$repro.repeats)
  
  rep.df = rbindlist(rep.df)
  rep.df = merge(rep.df,
                 rep.ests[, .(Effect.Index, RepSet, rep.sample.size,
                              Original.Effect.Size, CI.low, CI.high)],
                 by = "Effect.Index")
  
  return (rep.df)
}

# Goes over the replication data frame with evaluations and calculates overall summaries for
#   each of the measures present there.
reproducibility.rate = function (rates.df, n.sample = -1) {
  # Get a sample if required
  if (n.sample > 0) {
    eff.sample = sample(x = unique(rates.df$Effect.Index), size = n.sample, replace = F)
    target.df = rates.df[Effect.Index %in% eff.sample]
  }
  
  # Summarises
  # ...
}

# For each experiment, calculates reproducibility measures
calc.rep.measures = function(types = c("Orig-in-MA-PI", "MA-SSS", "VOTE-SSS-PAR", "VOTE-SSS-NPAR", "Orig-in-MA-CI", "MA-in-Orig-CI", "Orig-in-CI-3"), min.effect.of.interest) {
    
  # For each set of experiments, calculates each reproducibility measure in types
  rates.df = target.df[, evaluate.exp.rep(c(.BY, .SD), types = types,
                         min.effect.of.interest = min.effect.of.interest),
    by = .(Effect.Index, RepSet)]
  
  return (rates.df)
}

# Computes success or failure in a replication according to a give criterium
reproducibility.success = function (rep.exps, MA, type) {
  if (type == "Orig-in-MA-PI") {
    longtype = "original estimate is within the PI of the meta analysis"
    value = rep.exps$Original.Effect.Size[1] > MA$pred$cr.lb &
      rep.exps$Original.Effect.Size[1] < MA$pred$cr.ub
  }
  
  if (type == "Orig-in-MA-CI") {
    longtype = "original estimate is within the CI of the meta analysis"
    value = rep.exps$Original.Effect.Size[1] > MA$pred$ci.lb &
      rep.exps$Original.Effect.Size[1] < MA$pred$ci.ub
  }
  
  if (type == "MA-SSS") {
    longtype = "meta analysis is significant and in the same sense as the original"
    value = rep.exps$Original.Effect.Size[1] / MA$m$beta[[1]] > 0 & MA$m$pval < 0.05
  }
  
  if (type == "VOTE-SSS") {
    longtype = "majority of individual replications are significant and in the same sense"
    value = mean(rep.exps$p.value < 0.05) >= 0.5
  }
  
  if (type == "MA-in-Orig-CI") {
    longtype = "meta analysis point estimate is within the CI of the original"
    value = MA$beta[[1]] > rep.exps$CI.low[1] & MA$beta[[1]] < rep.exps$CI.high[1]
  }
  
  data.frame(
    Type = type,
    LongType = longtype,
    Success = value
  )
}

# Computes many types of reproducibility measures from a set of replications
evaluate.exp.rep = function (rep.exps, types, min.effect.of.interest) {
  
  MA = rep.exps[,run.ma(MeanControl, SDControl, Sample.Size,
                       MeanTreated, SDTreated, Sample.Size)]

  result = ldply(types, reproducibility.success(rep.exps, MA, type))
  original.estimate = rep.exps$Original.Effect.Size[1]
  real.effect = rep.exps$Real.Effect.Size[1]
  
  # By criteria
  # FP (reproduced but is not real)
  result$FP = result$Success & abs(original.estimate) < min.effect.of.interest
  # FN (didn't reproduce but is real)
  result$FN = !result$Success & abs(original.estimate) >= min.effect.of.interest
  # TP (reproduced and is real)
  result$TP = result$Success & abs(original.estimate) >= min.effect.of.interest
  # TN (didn't reproduce and is not real)
  result$TN = !result$Success & abs(original.estimate) < min.effect.of.interest
  
  # Cast to long format
  result = melt(result, id.vars = c("Type", "LongType", "Original.Effect.Size"))
  colnames(result) = c("Type", "LongType", "Original.Effect.Size", "Measure", "Value")
  
  # Exaggeration (MA x Original)
  result = rbind(
    result, data.frame(Type = NA, LongType = NA,
                       Measure = "Exaggeration (MA x Original)",
                       Value = ifelse(MA$m$beta[[1]] / original.estimate <= 0,
                                      NA, MA$m$beta[[1]] / original.estimate)))
  # Exaggeration (MA x Real)
  result = rbind(
    result, data.frame(Type = NA, LongType = NA,
                       Measure = "Exaggeration (MA x Real)",
                       Value = ifelse(MA$m$beta[[1]] / real.effect <= 0,
                                      NA, MA$m$beta[[1]] / real.effect)))
  # Signal (MA x Original)
  result = rbind(
    result, data.frame(Type = NA, LongType = NA,
                       Measure = "Signal Error (MA x Original)",
                       Value = MA$m$beta[[1]] / original.estimate <= 0))
  # Signal (MA x Real)
  result = rbind(
    result, data.frame(Type = NA, LongType = NA,
                       Measure = "Signal Error (MA x Real)",
                       Value = MA$m$beta[[1]] / real.effect <= 0))
  
  result$EffectIndex = rep.exps$EffectIndex[1]
  result$RepSet = rep.exps$RepSet[1]
  
  result
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
    # PI = list(PI.lower = pred$cr.lb, PI.upper = pred$cr.ub, beta = m$beta[[1]])
    # signif = list(eff = pred$pred,  se = pred$se, zval = m$zval, signif = abs(pred$pred) - 1.96 * pred$se > 0)
    # beta = list(beta = m$beta[[1]])
    return (pred, m)
  }, error = function(e) {
    return (NULL)
  })
}