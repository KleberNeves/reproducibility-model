library(data.table)
library(metafor)
library(tidyr)

# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.rep.evaluation.tests = function(min.effect, repro.detect) {
  
  rep.results.df = calc.rep.measures(min.effect.of.interest = min.effect, repro.detect = repro.detect)
  
  df = data.frame()
  
  df = rbind(df, reproducibility.rate(rep.results.df, n.sample = -1))
  df = rbind(df, reproducibility.rate(rep.results.df, n.sample = 20))
  
  df
}

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
    exps = map(1:separate.reps, function (rep_i, effect.index, rep.input) {
      exp = perform.experiment(effect.index, rep.input)
      exp$RepSet = rep_i
      exp
    }, effect.index, rep.input)
    
    rbindlist(exps)
  }

  # Runs a number of replications for each effect 
  rep.df = map(
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
reproducibility.rate = function (rep.results.df, n.sample = -1) {
  
  # Get a sample if required
  if (n.sample > 0) {
    eff.sample = sample(x = unique(rep.results.df$Effect.Index), size = n.sample, replace = F)
    sample.results.df = rep.results.df[Effect.Index %in% eff.sample]
  } else {
    sample.results.df = rep.results.df
  }
  
  # Summarises replication results (mean, min, max) for the sample
  rep.rates = sample.results.df %>% filter(!is.na(Type)) %>%
    pivot_wider(id_cols = c("Effect.Index", "Type", "LongType", "RepSet"),
                names_from = "Measure", values_from = "Value")
  walk(5:ncol(rep.rates), function(x) { rep.rates[[x]] <<- as.logical(rep.rates[[x]]) })
  
  # overall replication rate, PPV/precision, recall
  rep.rates = data.table(rep.rates)
  rep.rates = rep.rates[, .(
    ReproRate = mean(Success),
    Specificity_AM = sum(TN_AM) / (sum(TN_AM) + sum(FP_AM)),
    Sensitivity_AM = sum(TP_AM) / (sum(TP_AM) + sum(FN_AM)),
    Specificity_B = sum(TN_B) / (sum(TN_B) + sum(FP_B)),
    Sensitivity_B = sum(TP_B) / (sum(TP_B) + sum(FN_B)),
    Specificity_P = sum(TN_P) / (sum(TN_P) + sum(FP_P)),
    Sensitivity_P = sum(TP_P) / (sum(TP_P) + sum(FN_P))
  ), by = .(RepSet, Type, LongType)]
  
  # Prevalences in the sample
  cast.prev = function (x, measure) {
    x = x %>% filter(Measure == rlang::as_name(enquo(measure))) %>% select(-Type, -LongType) %>%
      pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                  names_from = "Measure", values_from = "Value") %>%
    mutate({{ measure }} := as.logical({{ measure }}))
    data.table(x)
  }
  
  temp.df = cast.prev(sample.results.df, Is.Above.Min)
  other.rates.am = temp.df[, .(Prev_Sample_AM = mean(Is.Above.Min)),
                           by = .(RepSet)]
  
  temp.df = cast.prev(sample.results.df, Is.Biased)
  other.rates.b = temp.df[, .(Prev_Sample_B = mean(Is.Biased)),
                          by = .(RepSet)]
  
  temp.df = cast.prev(sample.results.df, Is.Precise)
  other.rates.p = temp.df[, .(Prev_Sample_P = mean(Is.Precise)),
                          by = .(RepSet)]
  
  # Prevalence of the literature (whole literature, for reference)
  temp.df = cast.prev(rep.results.df, Is.Above.Min)
  other.rates.aml = temp.df[, .(Prev_Lit_AM = mean(Is.Above.Min)),
                                by = .(RepSet)]
  
  temp.df = cast.prev(rep.results.df, Is.Biased)
  other.rates.bl = temp.df[, .(Prev_Lit_B = mean(Is.Biased)),
                               by = .(RepSet)]
  
  temp.df = cast.prev(rep.results.df, Is.Precise)
  other.rates.pl = temp.df[, .(Prev_Lit_P = mean(Is.Precise)),
                               by = .(RepSet)]
  
  # General measures (criteria-independent), for the sample
  temp.df = sample.results.df %>% filter(is.na(Type) & !(Measure %in% c("Is.Above.Min", "Is.Biased", "Is.Precise"))) %>%
    select(-Type, -LongType) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value")
  walk(c(3,4,7,8), function(x) { temp.df[[x]] <<- as.numeric(temp.df[[x]]) })
  walk(c(5,6,9,10), function(x) { temp.df[[x]] <<- as.logical(temp.df[[x]]) })
  
  # Exaggeration and signal error rate
  temp.df = data.table(temp.df)
  error.rates = temp.df[,
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
  other.rates.am = pivot_longer(other.rates.am, cols = -"RepSet")
  other.rates.aml = pivot_longer(other.rates.aml, cols = -"RepSet")
  other.rates.b = pivot_longer(other.rates.b, cols = -"RepSet")
  other.rates.bl = pivot_longer(other.rates.bl, cols = -"RepSet")
  other.rates.p = pivot_longer(other.rates.p, cols = -"RepSet")
  other.rates.pl = pivot_longer(other.rates.pl, cols = -"RepSet")
  rep.rates = pivot_longer(rep.rates, cols = -c("RepSet", "Type", "LongType"))
  
  other.rates = rbind(error.rates, other.rates.am, other.rates.aml, other.rates.b, other.rates.bl, other.rates.p, other.rates.pl)
  other.rates$Type = NA
  other.rates$LongType = NA
  
  final.rates = rbind(rep.rates, other.rates)
  final.rates$N = ifelse(n.sample == -1, "All", as.character(n.sample))
  return (final.rates)
}

# For each experiment, calculates reproducibility measures
calc.rep.measures = function(types = c("Orig-in-RMA-PI", "RMA-SSS", "Orig-in-RMA-CI", "RMA-in-Orig-CI", "FMA-SSS", "Orig-in-FMA-CI", "FMA-in-Orig-CI", "VOTE-SSS"), min.effect.of.interest, repro.detect) {
  # For each set of experiments, calculates each reproducibility measure in types
  rates.df = replications.df[, evaluate.exp.rep(.SD, types = types,
                                                min.effect.of.interest = min.effect.of.interest, repro.detect = repro.detect),
                             by = .(Effect.Index, RepSet)]
  return (rates.df)
}

# Computes many types of reproducibility measures from a set of replications
evaluate.exp.rep = function (rep.exps, types, min.effect.of.interest, repro.detect) {
  
  RMA = with(rep.exps, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size, type = "RE"))
  
  FMA = with(rep.exps, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size, type = "FE"))
  
  result = map_dfr(types, reproducibility.success, rep.exps = rep.exps,
                 RMA = RMA, FMA = FMA)
  original.estimate = rep.exps$Original.Effect.Size[1]
  real.effect = rep.exps$Real.Effect.Size[1]
  reproduced = result$Success
  biased = rep.exps$Biased[1]
  is.above.min = abs(real.effect) >= min.effect.of.interest
  is.precise = abs(real.effect - original.estimate) <= repro.detect
  
  get_2x2_table = function (is_real, real_type) {
    df = data.frame(
      # False positives
      FP = reproduced & !is_real,
      # False negatives
      FN = !reproduced & is_real,
      # True positives
      TP = reproduced & is_real,
      # True negatives
      TN = !reproduced & !is_real
    )
    colnames(df) = paste(colnames(df), real_type, sep = "_")
    df
  }
  
  result = cbind(result,
                 get_2x2_table(is.above.min, "AM"), 
                 get_2x2_table(biased, "B"), 
                 get_2x2_table(is.precise, "P")) 
  
  # Cast to long format (seems convoluted, there should be a better way)
  result = pivot_longer(result, cols = -c("Type", "LongType"))
  colnames(result) = c("Type", "LongType", "Measure", "Value")
  result$Value = as.character(result$Value)
  
  result = rbind(result,
                 # Real Effect?
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Is.Above.Min",
                            Value = is.above.min),
                 
                 # Is Biased?
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Is.Biased",
                            Value = biased),
                 
                 # Well Estimated?
                 data.frame(Type = NA, LongType = NA,
                            Measure = "Is.Precise",
                            Value = is.precise),
                 
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
