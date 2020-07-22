# Returns a data frame with the results of replications
perform.replications = function(input, rep.power = -1, n.reps = 3) {
  # Filters the published estimates
  rep.ests = estimates.df[Published == T & p.value <= Alpha]
  
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
  replicate.exp = function(effect.index, rep.df) {
    rep.input = sanitize_shiny_input(input)
    rep.input$typical.sample.size = rep.ests[Effect.Index == effect.index, rep.sample.size]
    xp = perform.experiment(effect.index, rep.input, -1)
  }
  
  # Runs a number of replications for each effect
  rep.df = lapply(rep(rep.ests$Effect.Index, n.reps), replicate.exp, rep.df)
  rep.df = rbindlist(rep.df)
  rep.df = merge(rep.df, rep.ests[,.(Effect.Index, rep.sample.size, Original.Effect.Size, Original.p.value)], by = "Effect.Index")
  
  return (rep.df)
}

# Computes many types of reproducibility measures from a set of replications
reproducibility.rate = function(master.rep.df, type, input, n.sample = -1, measure.name = "", rep_i = 0) {
  
  # type is one of:
  # "BRI" (Brazilian Reproducibility Initiative)
  # "SSS" (significant, in the same sense)
  # "ST" (Small Telescopes)
  # n.sample is the number of experiments selected to evaluate (-1 means all experiments)
  
  n.experiments = ifelse(n.sample > 0, n.sample, nrow(master.rep.df)/3)
  
  # Pick the sample of experiments: sample is alawys picked in order of the effect index, to assure that every measure will be calculated on the same sample
  if (length(unique(master.rep.df$Effect.Index)) < n.experiments) {
    cat("Not enough experiments to reproduce!")
    return (-1)
  }
  
  # eff.sample = sort(unique(master.rep.df$Effect.Index))[1:n.experiments]
  eff.sample = sample(x = unique(master.rep.df$Effect.Index), size = n.experiments, replace = F)
  rep.df = master.rep.df[Effect.Index %in% eff.sample]
  rep.ests = estimates.df[Effect.Index %in% eff.sample]
  
  rep.ests = merge(rep.ests, rep.df %>% select("Effect.Index","rep.sample.size"), by = "Effect.Index")
  rep.ests = rep.ests[!duplicated(Effect.Index)]
  
  
  # Depending on the type of initiative, calculate the reproducibility rate
  
  # If BRI, do a meta analysis of the three replications and extract a prediction interval
  if (type == "BRI") {
    rep.df = rep.df[, run.ma(MeanControl, SDControl, Sample.Size,
                             MeanTreated, SDTreated, Sample.Size,
                             what = "PI"), Effect.Index]
    
    # rep.ests = rep.ests[!duplicated(Effect.Index)]
    
    rep.ests = merge(rep.ests, rep.df, by = "Effect.Index") %>%
      select(Effect.Index, Real.Effect.Size, Estimated.Effect.Size, PI.lower, PI.upper)
    rep.ests$Reproduced = rep.ests$Estimated.Effect.Size < rep.ests$PI.upper &
      rep.ests$Estimated.Effect.Size > rep.ests$PI.lower
    
    
    # If SSS, check if it is in the same sense of the original, and significant
  } else if (type == "SSS") {
    rep.res = rep.df[, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size,
                              what = "S"), Effect.Index]
    
    rep.ests = merge(rep.ests, rep.res, by = "Effect.Index") %>%
      select(Effect.Index, Real.Effect.Size, Estimated.Effect.Size, eff, signif)
    rep.ests$Reproduced = rep.ests$signif &
      rep.ests$eff / rep.ests$Estimated.Effect.Size > 0
    
    
    # If ST check if it is significantly different from a small effect (defined as an effect that the original study has 33% power to detect) in a one sided T-test at 5% level
    # The Small Telescopes approach recommends using a replication sample that is always 2.5 times the original sample size, as this gives about 80% power to  reject a "small effect"
  } else if (type == "ST") {
    
    #returns power calculation results for t-test at 33% power using original estimates 
    small_telescope_d33 = mapply(pwr.t.test, n = rep.ests$Sample.Size, power = 1/3)
    small_telescope_d33 = as.data.frame(small_telescope_d33)
    #transposes dataframe to correct rows vs columns
    small_telescope_d33 = as.data.frame(t(small_telescope_d33))
    #adds d33 to rep.ests
    rep.ests$small_telescope_d33 = t(as.data.frame(small_telescope_d33$d))
    
    rep.res = rep.df[, run.ma(MeanControl, SDControl, Sample.Size,
                              MeanTreated, SDTreated, Sample.Size,
                              what = "S"), Effect.Index]
    
    rep.ests = merge(rep.ests, rep.res, by = "Effect.Index") %>%
      select(Effect.Index, Real.Effect.Size, Estimated.Effect.Size, eff, se, rep.sample.size, signif, zval, small_telescope_d33)
    
    rep.ests$Reproduced = rep.ests$eff / rep.ests$Estimated.Effect.Size > 0 &
      pt(abs(rep.ests$zval), df = rep.ests$rep.sample.size - 1, ncp = rep.ests$small_telescope_d33 * sqrt(rep.ests$rep.sample.size / 2)) > 0.05
    rep.ests$small_telescope_d33 = NULL # GAMBIARRA, VER COM O PEDRO DEPOIS O QUE HOUVE ... PARECE QUE A COLUNA d_33 É UM DATA FRAME? ...
  }
  
  tried.true = rep.ests[abs(Real.Effect.Size) >= input$min.effect.of.interest, .N]
  tried.false = rep.ests[abs(Real.Effect.Size) < input$min.effect.of.interest, .N]
  FP = rep.ests[Reproduced == T & abs(Real.Effect.Size) < input$min.effect.of.interest, .N]
  FN = rep.ests[Reproduced == F & abs(Real.Effect.Size) >= input$min.effect.of.interest, .N]
  TP = rep.ests[Reproduced == T & abs(Real.Effect.Size) >= input$min.effect.of.interest, .N]
  TN = rep.ests[Reproduced == F & abs(Real.Effect.Size) < input$min.effect.of.interest, .N]
  
  reproduced = rep.ests[Reproduced == T, .N]
  
  rep.ests$Mfactor = ifelse(
    rep.ests$Real.Effect.Size / rep.ests$Estimated.Effect.Size > 0,
    rep.ests$Estimated.Effect.Size / rep.ests$Real.Effect.Size, NA)
  # filter for abs(Estimated.Effect.Size) >= Min.Interesting.Effect?
  
  # browser()
  
  return (
    data.frame(Measure = measure.name,
               Value = reproduced / n.experiments,
               TypeIErrorRate = FP / tried.false,
               TypeIIErrorRate = FN / tried.true,
               MFactor = mean(rep.ests$Mfactor, na.rm = T),
               TypeSErrorRate = mean(rep.ests$Real.Effect.Size / rep.ests$Estimated.Effect.Size < 0, na.rm = T),
               PPV = TP / (TP + FP),
               NPV = TN / (TN + FN),
               # Accuracy = ,
               `Published Effects Only?` = T,
               Sample.Size = n.sample,
               Rep = rep_i)
  )
}

# Runs a meta analysis given the means, SDs and Ns, then returns what is asked.
run.ma = function(mean_control, sd_control, n_control, mean_treated, sd_treated, n_treated, what) {
  
  ess = escalc(measure = "MD", m1i = as.numeric(mean_treated), 
               m2i = as.numeric(mean_control), sd1i = as.numeric(sd_treated), 
               sd2i = as.numeric(sd_control), n1i = as.numeric(n_treated), 
               n2i = as.numeric(n_control))
  tryCatch({
    m = rma(yi=yi, vi=vi, data = ess, measure = "MD", method = "REML", control=list(maxiter=1000))
    pred = predict.rma(m, level = 0.95, digits = 3)
    PI = list(PI.lower = pred$cr.lb, PI.upper = pred$cr.ub, beta = m$beta[[1]])
    signif = list(eff = pred$pred,  se = pred$se, zval = m$zval, signif = abs(pred$pred) - 1.96 * pred$se > 0)
    beta = list(beta = m$beta[[1]])
  }, error = function(e) {
    PI = list(PI.lower = -999, PI.upper = 999)
    signif = list(eff = 999, signif = FALSE)
    beta = list(beta = 999)
  })
  
  if (what == "PI") { PI } else if (what == "S") { signif } else if (what == "EST") { beta } else if (what == "RMA") { m }
}