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


# Helper functions for areas under curves
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
