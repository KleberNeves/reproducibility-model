library(data.table)
library(metafor)
library(tidyr)
library(BayesFactor)

# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.rep.evaluation.tests = function(input) {
  
  rep.results.df = calc.rep.measures(input)
  
  df = map_dfr(input$repro.sample, function (x) {
    reproducibility.rate(rep.results.df, n.sample = x)
  })
  
  df
}

# Returns a data frame with the results of replications
perform.replications = function(input) {
  # Filters the published estimates
  rep.ests = estimates.df[Published == T & p.value <= Alpha]
  
  # Number of replications per experiment
  n.reps = 3
  if (!is.null(input$repro.exps)) {
    n.reps = input$repro.exps
  }
  
  # Calculate N for the desired power for the original estimate (not the real effect)
  if (input$repro.power < 1) {
    calc.n = function (eff, sd, wanted.pwr, alpha) {
      pw = tryCatch({
        power.t.test(delta = eff, sd = sd, sig.level = alpha, power = wanted.pwr)$n
      }, error = function(e) {
        # Error means N < 2, i.e. very large effect compared to the SEM
        3
      })
      return (ceiling(pw))
    }
    
    rep.ests$rep.sample.size = mapply(calc.n,
                                      rep.ests$Estimated.Effect.Size,
                                      rep.ests$Estimated.Pooled.SD,
                                      input$repro.power, input$alpha.threshold)
  } else {
    rep.ests$rep.sample.size = input$typical.sample.size * input$repro.power
  }
  
  # Function to perform a single replication, for all the rep sets
  # Receives an effect index, and performs an experiment, adding it to the rep.df
  replicate.exp = function(effect.index, rep.df, separate.reps) {
    rep.input = sanitize_shiny_input(input)
    rep.input$typical.sample.size = rep.ests[Effect.Index == effect.index, rep.sample.size]
    
    # Runs separate repro repeats and saves that information so that we can estimate uncertainty later
    exps = map(1:separate.reps, function (rep_i, effect.index, rep.input) {
      exp = perform.experiment(effect.index, rep.input)
      exp$IsReplication = T
      exp$RepSet = rep_i
      exp
    }, effect.index, rep.input)
    
    rbindlist(exps)
  }
  
  # Runs a number of replications for each effect 
  rep.df = map(rep(rep.ests$Effect.Index, n.reps),
    replicate.exp, rep.df, separate.reps = input$repro.repeats
  )
  
  # To each rep set, add the original effect
  orig.df = map(1:input$repro.repeats, function (i) {
    exp = rep.ests[, !"rep.sample.size"]
    exp$IsReplication = F
    exp$RepSet = i
    exp
  })
  
  rep.df = rbindlist(c(rep.df, orig.df))
  
  return (rep.df)
}

# Goes over the replication data frame with evaluations and calculates overall summaries for each of the measures present there.
reproducibility.rate = function (rep.results.df, n.sample) {
  # browser()
  # Build the sample of replications/effects to be considered
  eff.sample = sample(x = unique(rep.results.df$Effect.Index),
                      size = n.sample, replace = F)
  
  sample.results.df = rep.results.df[Effect.Index %in% eff.sample]
  
  # Summarises replication results (mean, min, max) for the sample
  rep.rates = sample.results.df %>% filter(!is.na(Type)) %>%
    pivot_wider(id_cols = c("Effect.Index", "Type", "RepSet"),
                names_from = "Measure", values_from = "Value") %>%
    mutate(across(c(Success, TP, TN, FP, FN), as.logical))
  
  # overall replication rate, PPV/precision, recall
  rep.rates = data.table(rep.rates)
  rep.rates = rep.rates[, .(
    ReproRate = mean(Success),
    TN = sum(TN),
    TP = sum(TP),
    FP = sum(FP),
    FN = sum(FN)
  ), by = .(RepSet, Type)]
  
  cast.prev = function (x, measure) {
    x = x %>% filter(Measure == measure) %>% select(-Type) %>%
      pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                  names_from = "Measure", values_from = "Value")
    x[, measure] = as.logical(x[, measure, drop = T])
    data.table(x)
  }
  
  # Prevalences in the sample
  temp.df = cast.prev(sample.results.df, "Is.Real")
  other.rates.sample = temp.df[, .(Prev_Sample = mean(Is.Real, na.rm = T)), by = .(RepSet)]
  
  # Prevalence of the literature (whole literature, for reference)
  temp.df = cast.prev(rep.results.df, "Is.Real")
  other.rates.lit = temp.df[, .(Prev_Lit = mean(Is.Real, na.rm = T)), by = .(RepSet)]
  
  # Exaggeration and signal error rate, for the sample
  temp.df = sample.results.df %>% filter(is.na(Type) & Measure != "Is.Real") %>%
    select(-Type) %>%
    pivot_wider(id_cols = c("Effect.Index", "RepSet"),
                names_from = "Measure", values_from = "Value") %>%
    mutate(across(starts_with("Exaggeration"), as.numeric)) %>%
    mutate(across(starts_with("Signal"), as.logical)) %>%
    data.table()
  
  error.rates = temp.df[,
                        .(
                          Exaggeration_RMA_x_Original = median(`Exaggeration (RMA x Original)`,
                                                               na.rm = T),
                          Exaggeration_RMA_x_Real = median(`Exaggeration (RMA x Real)`,
                                                           na.rm = T),
                          Signal_Error_RMA_x_Original = mean(`Signal Error (RMA x Original)`,
                                                             na.rm = T),
                          Signal_Error_RMA_x_Real = mean(`Signal Error (RMA x Real)`,
                                                         na.rm = T)
                        ),
                        by = .(RepSet)]
  
  # Make long, bind and return
  error.rates = pivot_longer(error.rates, cols = -"RepSet")
  other.rates.sample = pivot_longer(other.rates.sample, cols = -"RepSet")
  other.rates.lit = pivot_longer(other.rates.lit, cols = -"RepSet")
  rep.rates = pivot_longer(rep.rates, cols = -c("RepSet", "Type"))
  
  other.rates = rbind(error.rates, other.rates.sample, other.rates.lit)
  other.rates$Type = NA
  
  final.rates = rbind(rep.rates, other.rates)
  final.rates$N = n.sample
  final.rates$Nprop = round(n.sample / length(unique(rep.results.df$Effect.Index)), 2)
  return (final.rates)
}

# For each experiment, calculates all the reproducibility measures
calc.rep.measures = function(input) {
  # For each set of experiments, calculates each reproducibility measure in types
  rates.df = replications.df[, evaluate.exp.rep(.SD, input = input),
                             by = .(Effect.Index, RepSet)]
  return (rates.df)
}

# Computes many types of reproducibility measures from a set of replications
evaluate.exp.rep = function (rep.exps, input) {
  
  if (nrow(rep.exps[IsReplication == T,]) > 1) {  
    RMA = with(rep.exps[IsReplication == T,],
               run.ma(MeanControl, SDControl, Sample.Size,
                      MeanTreated, SDTreated, Sample.Size, type = "RE"))
    
    FMA = with(rep.exps[IsReplication == T,],
               run.ma(MeanControl, SDControl, Sample.Size,
                      MeanTreated, SDTreated, Sample.Size, type = "FE"))
    
    rep_estimate = RMA$m$beta[[1]]; rep_p = RMA$m$pval
  } else {
    RMA = NA; FMA = NA
    rep_estimate = rep.exps[IsReplication == T, Estimated.Effect.Size][1]
    rep_p = rep.exps[IsReplication == T, p.value][1]
  }
  
  CMA = with(rep.exps,
             run.ma(MeanControl, SDControl, Sample.Size,
                    MeanTreated, SDTreated, Sample.Size, type = "RE"))
  
  result = reproducibility.success(rep.exps, RMA, FMA, CMA)
  
  original.estimate = rep.exps[IsReplication == F, Estimated.Effect.Size][1]
  real.effect = rep.exps[IsReplication == F, Real.Effect.Size][1]
  reproduced = result$Success
  
  if (input$fixed.prev.mode == "none" | input$fixed.prev.mode == "above minimum of interest") {
    is_real = abs(real.effect) >= input$min.effect.of.interest
  } else if (input$fixed.prev.mode == "non-biased") {
    is_real = !rep.exps$Biased[1]
  } else if (input$fixed.prev.mode == "precision") {
    is_real = abs(real.effect - original.estimate) <= input$repro.detect
  }
  
  get_2x2_table = function (is_real) {
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
    df
  }
  
  result = cbind(result, get_2x2_table(is_real)) 
  
  # Cast to long format (seems convoluted, there should be a better way)
  result = pivot_longer(result, cols = -c("Type"), names_to = "Measure", values_to = "Value")
  result$Value = as.character(result$Value)
  
  result = rbind(result,
                 # Real Effect?
                 data.frame(Type = NA,
                            Measure = "Is.Real",
                            Value = is_real),
                 
                 # Exaggeration (RMA x Original)
                 data.frame(Type = NA,
                            Measure = "Exaggeration (RMA x Original)",
                            Value = ifelse(rep_estimate / original.estimate <= 0,
                                           NA, rep_estimate / original.estimate)),
                 # Exaggeration (RMA x Real)
                 data.frame(Type = NA,
                            Measure = "Exaggeration (RMA x Real)",
                            Value = ifelse(rep_estimate / real.effect <= 0,
                                           NA, rep_estimate / real.effect)),
                 # Signal (RMA x Original)
                 data.frame(Type = NA,
                            Measure = "Signal Error (RMA x Original)",
                            Value = rep_estimate / original.estimate <= 0 &
                              rep_p < 0.05),
                 # Signal (RMA x Real)
                 data.frame(Type = NA,
                            Measure = "Signal Error (RMA x Real)",
                            Value = rep_estimate / real.effect <= 0 &
                              rep_p < 0.05)
  )
  
  result
}

# Computes success or failure in a replication according to a given criterion
reproducibility.success = function (comb.exps, RMA, FMA, CMA) {
  
  rep.exps = comb.exps[IsReplication == T, ]
  orig.exp = comb.exps[IsReplication == F, ]
  
  d33 = power.t.test(
    n = orig.exp$Sample.Size, delta = NULL,
    sd = 1, sig.level = orig.exp$Alpha,
    power = 0.33
  )$delta
  
  sdps = sqrt((rep.exps$SDControl ^ 2 + rep.exps$SDTreated ^ 2) / 2)
  t_scores = rep.exps$Estimated.Effect.Size / sdps / sqrt(2 / rep.exps$Sample.Size)
  bf = (meta.ttestBF(t = t_scores, n1 = rep.exps$Sample.Size, n2 = rep.exps$Sample.Size, rscale = 1))@bayesFactor$bf
  
  types = c("VOTE_SSS_005", "VOTE_SSS_0005", "FMA_SSS_005", "FMA_SSS_0005", "RMA_SSS_005", "RMA_SSS_0005", "ORIG_IN_RMA_PI", "ORIG_IN_FMA_CI", "REP_IN_ORIG_CI", "CMA_SSS_005", "CMA_SSS_0005", "SMALL_TELESCOPE", "BF_3", "BF_10")
  
  successes = c(
    # Simple majority voting by significance (p < 0.05) and same sense
    mean(rep.exps$p.value < 0.05) >= 0.5,
    
    # Simple majority voting by significance (p < 0.005) and same sense
    mean(rep.exps$p.value < 0.005) >= 0.5,
    
    # Significance (p < 0.05) and same sense of a fixed-effects meta-analysis
    ifelse(
      any(is.na(FMA)), NA,
      orig.exp$Estimated.Effect.Size[1] / FMA$m$beta[[1]] > 0 & FMA$m$pval < 0.05
    ),
    
    # Significance (p < 0.005) and same sense of the fixed-effects meta-analysis
    ifelse(
      any(is.na(FMA)), NA,
      orig.exp$Estimated.Effect.Size[1] / FMA$m$beta[[1]] > 0 & FMA$m$pval < 0.005
    ),
    
    # Significance (p < 0.05) and same sense of a random-effects meta-analysis
    ifelse(
      any(is.na(RMA)), NA,
      orig.exp$Estimated.Effect.Size[1] / RMA$m$beta[[1]] > 0 & RMA$m$pval < 0.05
    ),
    
    # Significance (p < 0.005) and same sense of the random-effects meta-analysis
    ifelse(
      any(is.na(RMA)), NA,
      orig.exp$Estimated.Effect.Size[1] / RMA$m$beta[[1]] > 0 & RMA$m$pval < 0.005
    ),
    
    # Original estimate is within the prediction interval of the random-effects meta-analysis of the replications
    ifelse(
      any(is.na(RMA)), NA,
      orig.exp$Estimated.Effect.Size[1] > RMA$pred$cr.lb &
        orig.exp$Estimated.Effect.Size[1] < RMA$pred$cr.ub
    ),
    
    # Original estimate is within the confidence interval of the fixed-effects meta-analysis of the replications
    ifelse(
      any(is.na(FMA)), NA,
    orig.exp$Estimated.Effect.Size[1] > FMA$pred$ci.lb &
      orig.exp$Estimated.Effect.Size[1] < FMA$pred$ci.ub
    ),
    
    # Summary of the replications is within the confidence interval of the original estimate
    ifelse(
      any(is.na(RMA)),
      rep.exps$Estimated.Effect.Size[1] > orig.exp$CI.low[1] & rep.exps$Estimated.Effect.Size[1] < orig.exp$CI.high[1],
      RMA$m$beta[[1]] > orig.exp$CI.low[1] & RMA$m$beta[[1]] < orig.exp$CI.high[1]
    ),
    
    # Combined meta-analysis is significant (p < 0.05)
    orig.exp$Estimated.Effect.Size[1] / CMA$m$beta[[1]] > 0 & CMA$m$pval < 0.05,
    
    # Combined meta-analysis is significant (p < 0.005)
    orig.exp$Estimated.Effect.Size[1] / CMA$m$beta[[1]] > 0 & CMA$m$pval < 0.005,
    
    # Small telescopes
    ifelse(
      any(is.na(RMA)),
      rep.exps$Estimated.Effect.Size[1] > d33,
      RMA$m$beta[[1]] > d33
    ),
    
    # Bayes factor for the alternative against the null hypothesis is larger than 3
    bf >= 3,
    
    # Bayes factor for the alternative against the null hypothesis is larger than 10
    bf >= 10
  )
  
  success_df = tibble(Type = types, Success = successes)
  
  success_df
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
