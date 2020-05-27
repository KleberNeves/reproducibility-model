# Sets up the data.table objects that will hold the experimental results of the simulation. These are preallocated in memory for performance.
# Most of the model's work is on these tables
setup.model = function(input) {
  m = input$sim.end.value
  
  effects.df <<- data.table()
  
  estimates.df <<- data.table(
    Effect.Index = numeric(m),
    Real.Effect.Size = numeric(m),
    Estimated.Effect.Size = numeric(m),
    Estimated.Pooled.SD = numeric(m),
    MeanControl = numeric(m),
    SDControl = numeric(m),
    MeanTreated = numeric(m),
    SDTreated = numeric(m),
    p.value = numeric(m),
    CI.low = numeric(m),
    CI.high = numeric(m),
    Sample.Size = numeric(m),
    Power = numeric(m),
    Interlab.Variation = numeric(m),
    Measurement.Error = numeric(m),
    Alpha = numeric(m),
    How.It.Ends = character(m),
    End.Value = numeric(m),
    Replication.Chance = numeric(m),
    Positive.Bias = numeric(m),
    Min.Interesting.Effect = numeric(m),
    Scientist.Index = numeric(m),
    Biased = logical(m),
    Published = logical(m),
    Is.Replication = logical(m)
  )
  
  estimates.rowcount <<- 0
  
  eval.df <<- data.frame()
  
  # Generating a sample from the underlying distribution
  dist.n = 100000
  a.sample = abs(sample.from.dist(
    input$sdA, input$weightB, input$meanB, input$sdB, dist.n))
  
  # Calculated parameters
  # % Above minimum
  p2 = mean(a.sample >= input$min.effect.of.interest)
  
  # browser()
  # Calculate sample size in the input
  # In this version, you don't set the sample size, only the power
  # the other one will be calculated from the other
  
  # Mean effect size for the effects above the minimum of interest
  mes = mean(a.sample[a.sample >= input$min.effect.of.interest])
  if (is.nan(mes)) { mes = input$min.effect.of.interest }
  
  input$typical.sample.size = ceiling(power.t.test(
    n = NULL, delta = mes,
    sd = 1, sig.level = input$alpha.threshold,
    power = input$typical.power
  )$n)

  x = input
  
  print(input$typical.sample.size)
  print(input$typical.power)
  print(mes)
  x["loadDataFile"] = NULL
  x = as.data.frame(as.matrix(unlist(x)), stringsAsFactors = F) %>%
    tibble::rownames_to_column()
  colnames(x) = c("Parameter", "Value")
  
  x = rbind(x, data.frame(
    data.frame(Parameter = "% Above minimum", Value = signif(p2,2))
  ))
  
  param.df <<- x
  
  input
}

# Obtains a sample from a mixture distribution of two normals specified by 4 parameters
# normal A is centered at 0
# sdA: SD of normal A
# weightB: percentage of mixture that comes from normal B
# meanB, sdB: center and SD of normal B
# k: number of samples to obtain
sample.from.dist = function (sdA, weightB, meanB, sdB, k = 1) {
  smax = max(1000, k)
  dA = rnorm(ceiling(smax * (1 - weightB)), 0, sdA)
  dB = rnorm(ceiling(smax * weightB), meanB, sdB)
  # browser()
  d = sample(c(dA,dB), size = k)
  d
}

sample.dichotomous = function (sdAB, weightB = 0.5, k = 1) {
  sample.from.dist(
    sdA = sdAB, weightB = weightB, meanB = 1, sdB = sdAB, k = k)
}

sample.single = function (sdA, k = 1) {
  sample.from.dist(
    sdA = sdA, weightB = 0, meanB = 0, sdB = 0, k = k)
}

sample.double = function (sdB, weightB = 0.5, k = 1) {
  sample.from.dist(
    sdA = 0.1, weightB = weightB, meanB = 0, sdB = sdB, k = k)
}

# Generates a new effect to be investigated
generate.effect.size = function(p_sdA, p_weightB, p_meanB, p_sdB, eff.dist) {
  eff = sample.from.dist(p_sdA, p_weightB, p_meanB, p_sdB)
  
  ind = nrow(effects.df) + 1
  effects.df <<- rbindlist(list(
    effects.df, data.frame(Effect.Index = ind, Effect.Size = eff)
  ))
  
  return (ind)
}

# Picks from published effects table and finds the index of the effect in the real table
# Assumes the same lab in another time is a different lab (i.e. would have a different interlab variation)
pick.published.effect = function() {
  ind = sample(unique((estimates.df %>% filter(Published))$Effect.Index), 1)
  return (ind)
}

# Add interlab variation (picked from a normal distribution centered on the real effect)
add.interlab.variation = function(effect.index, interlab.variation) {
  real.eff = effects.df[[effect.index,"Effect.Size"]]
  eff = rnorm(n = 1, mean = real.eff, sd = interlab.variation)
  return (eff)
}

# Pick a sample from two normal distributions (control group has a mean of 0, exp group has a mean of "effect size post interlab variation")
# Returns a data.table with two columns (value, group)
get.samples = function(n, effect.size) {
  
  control.values = rnorm(n, 0, 1)
  exp.values = rnorm(n, effect.size, 1)
  
  dcontrol = data.frame(Group = "Control", Value = control.values)
  dexp = data.frame(Group = "Exp", Value = exp.values)
  
  d = rbindlist(list(dcontrol, dexp))
  return (d)
}

# Add measurement error on each element of the sample
get.measurements = function(exp.data, error) {
  measured = exp.data$Value + rnorm(nrow(exp.data), 0, error)
  return (measured)
}

# Performs the statistical test
statistical.test = function(exp.data, alpha) {
  # Perform a t test for the two groups
  m = t.test(x = exp.data[Group == "Exp", MeasuredValue],
             y = exp.data[Group == "Control", MeasuredValue],
             alternative = "two.sided", var.equal = T)
  # Returns a clean summary (broomed + whatever else I care about, e.g. effect size)
  m = tidy(m)
  
  m$significant = m$p.value <= alpha
  m$sd1 = exp.data[Group == "Exp", sd(MeasuredValue)]
  m$sd2 = exp.data[Group == "Control", sd(MeasuredValue)]
  n = nrow(exp.data) / 2
  m$pooled.sd = ((m$sd1 ^ 2 + m$sd2 ^ 2) / 2) ^ 0.5
  return (m)
}

# Returns whether or not to pick a new or published effect, depending on the parameter for replication incentive
is.replication = function(rep.chance) {
  return (runif(1,0,1) < rep.chance)
}

# Returns whether or not it is published, based on result and publication bias parameters
published = function(p.value, Alpha, negative.bias) {
  return (p.value <= Alpha | runif(1,0,1) < negative.bias)
}

# Check whether the condition to end the simulation was reached (sample size or number of effects)
reached.sim.end = function(max.value, criteria, alpha, perc = 1) {
  return (sim.end.tracking(criteria, alpha) >= max.value * perc)
}

# Tracks the value that determines the end of the simulation, depending on the criteria chosen by the user
sim.end.tracking = function(criteria, alpha) {
  if (criteria == "At a given number of published effects") {
    v = estimates.df[Published == T, .N]
  } else if (criteria == "At a given number of positive published effects") {
    v = estimates.df[Published == T & p.value <= alpha, .N]
  } else if (criteria == "At a given number of unique positive published effects") {
    v = length(unique(estimates.df[Published == T & p.value <= alpha, Effect.Index]))
  } else if (criteria == "At a given total sample size") {
    v = sum(estimates.df$Sample.Size)
  }
  v
}

# Function to simulate the whole experimental procedure, from data generation to analysis and publication
# Calls all the others functions; works as an outline of the main process in the model
perform.experiment = function(effect.index, input, scientist.id) {
  
  exp.effect = effects.df[Effect.Index == effect.index,]
  
  # Necessário ser uma coluna???
  exp.effect$effect.size.here = add.interlab.variation(effect.index, input$interlab.var)
  
  exp.data = get.samples(input$typical.sample.size, exp.effect$effect.size.here)
  
  # if (input$measure.error > 0) {
  #   exp.data$MeasuredValue = get.measurements(exp.data, input$measure.error)
  # } else {
    exp.data$MeasuredValue = exp.data$Value
  # }
  
  test = statistical.test(exp.data, input$alpha.threshold)
  
  new.estimate = data.table(
    Effect.Index = exp.effect$Effect.Index,
    Real.Effect.Size = exp.effect$Effect.Size,
    Estimated.Effect.Size = test$estimate1 - test$estimate2,
    Estimated.Pooled.SD = test$pooled.sd,
    MeanControl = test$estimate2,
    SDControl = test$sd2,
    MeanTreated = test$estimate1,
    SDTreated = test$sd1,
    p.value = test$p.value,
    CI.low = test$conf.low,
    CI.high = test$conf.high,
    
    # The list below are fixed parameters - I could reduce RAM and storage
    # use if this is not used later (replace columns with the input value)
    # Would need to save/load the input as well, when saving results, when
    # opening later for analysis. But doing it like this is silly.
    # Changing it would also make it simpler when adding new parameters
    Sample.Size = input$typical.sample.size,
    Power = input$typical.power,
    Interlab.Variation = input$interlab.var,
    Measurement.Error = 0,#input$measure.error,
    Alpha = input$alpha.threshold,
    How.It.Ends = input$how.sim.ends,
    End.Value = input$sim.end.value,
    Replication.Chance = input$rep.incentive,
    Positive.Bias = input$neg.incentive,
    Min.Interesting.Effect = input$min.effect.of.interest,
    
    Scientist.Index = scientist.id,
    Biased = F,
    Published = published(test$p.value, input$alpha.threshold, input$neg.incentive)
  )
  
  return (new.estimate)
}

# Executes the simulated scientist behavior
scientist.action = function(input, scientist.id) {
  # If simulation ended, return
  if (reached.sim.end(input$sim.end.value, input$how.sim.ends, input$alpha.threshold)) {
    return (1)
  }
  
  # Define whether it will be a replication or an original effect
  if (nrow(estimates.df %>% filter(Published)) > 0 & is.replication(input$rep.incentive)) {
    effect.index = pick.published.effect()
    xp.Is.Replication = T
  } else {
    effect.index = generate.effect.size(input$sdA, input$weightB, input$meanB, input$sdB)
    xp.Is.Replication = F
  }
  
  # Performs the experiment
  xp = perform.experiment(effect.index, input, scientist.id)
  
  # Bias
  ### With a given probability, if the result is not significant,
  ### change the estimated effect size to be just over the threshold
  
  if (xp$p.value > xp$Alpha & runif(1,0,1) < input$bias.level) {
    # Calc pooled SD, then make the treated mean be control mean + 2 SEM (using pooled SD)
    sem = xp$Estimated.Pooled.SD #/ sqrt(input$typical.sample.size)
    xp$Estimated.Effect.Size = qt(xp$Alpha / 2, df = 2 * input$typical.sample.size - 2) * sem * sign(xp$Estimated.Effect.Size)
    xp$MeanTreated = xp$MeanControl + xp$Estimated.Effect.Size
    xp$p.value = xp$Alpha - 0.01
    xp$Published = T
    xp$Biased = T
  }
  
  # Increases the row count, allocates memory if the original allocation is full
  # This is for performance reasons, to avoid growing the data frame (see R Circles of Hell)
  estimates.rowcount <<- estimates.rowcount + 1
  if (estimates.rowcount > nrow(estimates.df)) {
    tmp = estimates.df
    tmp[1:nrow(tmp),] = NA
    estimates.df <<- rbindlist(list(estimates.df, tmp), use.names = T, fill = T)
  }
  
  # Deu erro na linha abaixo
  # [1] "Generating the literature ..."
  # Error in set(estimates.df, i = as.integer(estimates.rowcount), j, xp[[j]]) : 
  #   Supplied 7 items to be assigned to 1 items of column 'Effect.Index'. The RHS length must either be 1 (single values are ok) or match the LHS length exactly. If you wish to 'recycle' the RHS please use rep() explicitly to make this intent clear to readers of your code.
  # In addition: There were 50 or more warnings (use warnings() to see the first 50)
  # browser()
  # Adds new row using data.table::set  
  for (j in 1:ncol(xp)) {
    set(estimates.df, i = as.integer(estimates.rowcount), j, xp[[j]])
  }
  # browser()
  
  # Esse experimento e todos os outros do mesmo Effect Index agora são replicações, caso esse seja
  estimates.df[Effect.Index == effect.index, Is.Replication := xp.Is.Replication]

  return (0)
}

# Runs actions for each scientist (might be obsolete/unnecessary)
run.iteration = function(input) {
  for (i in 1:input$n.scientists) {
    x = scientist.action(input, i)
    if (x != 0) {
      return (x)
    }
  }
  
  return (0)
}

# Runs a replicate of the simulation (nice to abstract the run.simulation function, however, it is a legacy from when we ran multiple simulations at once from the GUI)
run.replicate = function(input) {
  t = glue("Running the simulation ...")
  print(t); if (shiny_running) { showNotification(t, type = "message", duration = 8) }
  
  t = "Generating the literature ..."
  print(t); if (shiny_running) { showNotification(t, type = "default") }
  
  it = 0
  r25 = F; r50 = F; r75 = F
  while (it != 1) {
    it = run.iteration(input)
    
    # Progress feedback
    if (!r25 &
        reached.sim.end(input$sim.end.value, input$how.sim.ends, input$alpha.threshold, 0.25)) {
      t = "25% done ..."
      print(t); if (shiny_running) { showNotification(t, type = "default") }
      r25 = T
    }
    if (!r50 &
        reached.sim.end(input$sim.end.value, input$how.sim.ends, input$alpha.threshold, 0.5)) {
      t = "Halfway there ..."
      print(t); if (shiny_running) { showNotification(t, type = "default") }
      r50 = T
    }
    if (!r75 &
        reached.sim.end(input$sim.end.value, input$how.sim.ends, input$alpha.threshold, 0.75)) {
      t = "Almost done ..."
      print(t); if (shiny_running) { showNotification(t, type = "default") }
      r75 = T
    }

  }
}

# Main function, called from the interface to set the whole model running
run.simulation = function(input) {
  if (shiny_running) {
    Linput = reactiveValuesToList(input)
  } else {
    Linput = input
  }
  
  Linput$n.scientists = 1
  
  evdf = data.frame()
  
  Linput = setup.model(Linput)
  
  run.replicate(Linput)
  
  # summary.estimates.df <<- estimates.df %>% mutate(Significant = p.value < Alpha, True.Effect = abs(Real.Effect.Size) >= Min.Interesting.Effect) %>% group_by(Published, Significant, True.Effect) %>% summarise(N = n()) %>% mutate(Correct = (Significant == True.Effect)) %>% select(2,3,1,5,4)
  
  # Save the presynthesis results in another dataframe before synthesizing
  presynthesis.df <<- estimates.df[, .SD[1:estimates.rowcount]]
  estimates.df <<- synthesize(presynthesis.df)
  
  t = "Evaluating the literature ..."
  print(t); if (shiny_running) { showNotification(t, type = "default") }
  
  new.measures = make.evaluation.tests(Linput)
  
  if (input$scenarioName == "") {
    new.measures$scenario = Linput$scenario 
  } else {
    new.measures$scenario = Linput$scenarioName
  }
  
  evdf = rbind(evdf, new.measures)
  
  eval.df <<- evdf

  t = "... and ... Finished!"
  print(t); if (shiny_running) { showNotification(t, type = "error") }
}

# Synthesizes the literature, performing meta analysis of replications of the same effect when appropriate
synthesize = function(raw.df) {
  
  make.one.synthesis = function (ADT) {
    m = run.ma(ADT$MeanControl, ADT$SDControl, ADT$Sample.Size, ADT$MeanTreated, ADT$SDTreated, ADT$Sample.Size, what = "RMA")
    
    list(
      Estimated.Effect.Size = m$beta[[1]],
      p.value = m$pval,
      Published = T,#m$pval < ADT$Alpha,
      Biased = F,
      Is.Replication = F,
      Sample.Size = sum(ADT$Sample.Size),
      N.Studies = nrow(ADT)
    )
  }
  
  raw.df$N.Studies = 1L
  syn.df = raw.df[Is.Replication == T,
                  c("Estimated.Effect.Size","p.value","Published","Biased","Is.Replication","Sample.Size","N.Studies") := make.one.synthesis(.SD),
                  by = Effect.Index]
  
  syn.df = syn.df[!duplicated(Effect.Index)]
  syn.df = syn.df %>% select(-MeanControl, -SDControl, -MeanTreated, -SDTreated)
  
  return (syn.df)
}

# Measures to evaluate the corpus of estimates of effect sizes
# General function that calculates all rates, in all cases
make.evaluation.tests = function(input) {
  df = data.frame()
  
  t = "Error rates and discoveries ..."
  print(t); if (shiny_running) { showNotification(t, type = "default") }
  
  pb = T
  df = rbind(df, data.frame(Measure = "True Positives", Value = minimum.effect.count(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "True Positives", Value = minimum.effect.count(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "True Negatives", Value = true.negatives(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "True Negatives", Value = true.negatives(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df,data.frame(Measure = "False Positives", Value = typeI.error.rate(pb),
                           `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "False Positives", Value = typeI.error.rate(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "False Negatives", Value = typeII.error.rate(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "False Negatives", Value = typeII.error.rate(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Signal Errors", Value = typeS.error.rate(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Signal Errors", Value = typeS.error.rate(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Exaggeration Factor", Value = exaggeration.factor(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Exaggeration Factor", Value = exaggeration.factor(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Signal Errors (Effects > Min)", Value = typeS.error.rate.above.min(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Signal Errors (Effects > Min)", Value = typeS.error.rate.above.min(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Exaggeration Factor (Effects > Min)", Value = exaggeration.factor.above.min(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Exaggeration Factor (Effects > Min)", Value = exaggeration.factor.above.min(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Positive Predictive Value", Value = pos.pred.value(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Positive Predictive Value", Value = pos.pred.value(pb),
                            `Published Effects Only?` = pb))
  
  pb = T
  df = rbind(df, data.frame(Measure = "Negative Predictive Value", Value = neg.pred.value(pb),
                            `Published Effects Only?` = pb))
  pb = F
  df = rbind(df, data.frame(Measure = "Negative Predictive Value", Value = neg.pred.value(pb),
                            `Published Effects Only?` = pb))

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

# TYPE I/FALSE POSITIVE ERROR RATE (CONDITIONAL ON SIGNAL BEING RIGHT)
typeI.error.rate.signal = function(published.only = T) {
  if (published.only) {
    ests = estimates.df %>% filter(Published == T)
  } else {
    ests = estimates.df
  }
  
  if (nrow(ests) == 0) {
    return (NA)
  }
  
  # Wrong Signal
  ests$signalError = (ests$Real.Effect.Size / ests$Estimated.Effect.Size < 0)
  
  ests = ests %>% filter(!signalError)
  
  if (nrow(ests) == 0) {
    return (NA)
  }
  
  # False Positives
  ests$isError = (abs(ests$Real.Effect.Size) < ests$Min.Interesting.Effect) & (ests$p.value <= ests$Alpha)
  
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
  ests = estimates.df %>% filter(p.value < Alpha, abs(Estimated.Effect.Size) >= Min.Interesting.Effect)
  
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
### EVALUATION FUNCTIONS WITHIN THE BLOCK ABOVE


# Returns a data frame with the results of replications
perform.replications = function(input, rep.power, n.reps = 3, positive.only = T) {
  # rep.power = -1 means "use the same sample size as the original"
  
  # Filters the published estimates
  rep.ests = estimates.df[Published == T]
  
  if (positive.only) {
    rep.ests = rep.ests[p.value < Alpha]
  }
  
  rep.ests = rep.ests[!duplicated(Effect.Index)]
  
  # Calculate N for the desired power according to the original estimate (not the real effect)
  if (rep.power > 0) {
    calc.n = function (eff, sd, wanted.pwr, alpha) {
      pw = tryCatch(
        { power.t.test(delta = eff, sd = sd, sig.level = alpha, power = wanted.pwr)$n },
        error = function(e) { 2 # Error means N < 2, i.e. very large effect compared to the SE
        })
      return (ceiling(pw))
    }
    rep.ests$rep.sample.size = mapply(calc.n, rep.ests$Estimated.Effect.Size, rep.ests$Estimated.Pooled.SD, rep.power, input$alpha.threshold)
  } else {
    rep.ests$rep.sample.size = input$typical.sample.size
  }
  
  # Receives an effect index, and performs an experiment, adding it to the rep.df
  replicate.exp = function(effect.index, rep.df) {
    if (shiny_running) {
      rep.input = reactiveValuesToList(input)
    } else {
      rep.input = input
    }
    rep.input$typical.sample.size = rep.ests[Effect.Index == effect.index, rep.sample.size]
    xp = perform.experiment(effect.index, rep.input, -1)
  }
  
  rep.ests$Original.Effect.Size = rep.ests$Estimated.Effect.Size
  rep.ests$Original.p.value = rep.ests$p.value
  
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

area.under.tails = function(pfunction, params, threshold) {
  x = 2 * (
    do.call(pfunction, c(q = threshold, params)) -
    do.call(pfunction, c(q = 0, params))
  )
  return (1 - x)
}

area.under.normal.tails = function(spread, threshold) {
  2 * (1- pnorm(mean = 0, sd = spread, q = threshold))
}

# Helper to find the prevalence for a given distribution (tests many different parameters and returns an estimate of the prevalence of above minimum effects)
get.above.minimum.weight = function(weights, mint, sdA, meanB, sdB, dist.n = 20000) {
  m = abs(sapply(weights, sample.from.dist, sdA = sdA, sdB = sdB, meanB = meanB, k = dist.n)) >= mint
  m = as.list(colSums(m) / dist.n)
  names(m) = weights
  m
}

get.above.minimum.sd = function(sdAs, mint, dist.n = 20000) {
  m = abs(sapply(sdAs, sample.from.dist, weightB = 0, sdB = 0, meanB = 0, k = dist.n)) >= mint
  m = as.list(colSums(m) / dist.n)
  names(m) = sdAs
  m
}
