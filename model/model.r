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
    Alpha = numeric(m),
    Min.Interesting.Effect = numeric(m),
    Biased = logical(m),
    Published = logical(m)
  )
  
  estimates.rowcount <<- 0
  
  eval.df <<- data.frame()
  rep.eval.df <<- data.frame()
  
  # Generating a sample from the underlying distribution
  dist.n = 100000
  a.sample = abs(sample.from.dist(
    input$sdA, input$weightB, input$meanB, input$sdB, dist.n))
  
  # Calculated parameters
  # % Above minimum
  p2 = mean(a.sample >= input$min.effect.of.interest)
  
  # Calculate sample size in the input
  if (input$power.calculation == "mean of the true effects") {
    es_to_detect = mean(a.sample[a.sample >= input$min.effect.of.interest])
    if (is.nan(es_to_detect)) { es_to_detect = input$min.effect.of.interest }
  } else if (input$power.calculation == "minimum of interest") {
    es_to_detect = input$min.effect.of.interest
  } else if (input$power.calculation == "twice the minimum of interest") {
    es_to_detect = 2 * input$min.effect.of.interest
  } else {
    stop("Invalid way of calculating power for experiments.")
  }
  # Mean effect size for the effects above the minimum of interest
  input$typical.sample.size = ceiling(power.t.test(
    n = NULL, delta = es_to_detect,
    sd = 1, sig.level = input$alpha.threshold,
    power = input$typical.power
  )$n)
  
  x = input
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
generate.effect.size = function(p_sdA, p_weightB, p_meanB, p_sdB, interlab.variation = 0) {
  eff = sample.from.dist(p_sdA, p_weightB, p_meanB, p_sdB)
  ind = nrow(effects.df) + 1
  effects.df <<- rbindlist(list(
    effects.df, data.frame(Effect.Index = ind, Effect.Size = eff,
                  Effect.Variation = get.interlab.variation(interlab.variation))
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
get.interlab.variation = function(interlab.variation) {
   return (interlab.variation) #runif(n = 1, min = 0, max = interlab.variation))
}

# Add interlab variation (picked from a normal distribution centered on the real effect)
add.interlab.variation = function(effect.index, interlab.variation) {
  real.eff = effects.df[[effect.index,"Effect.Size"]]
  effvar = effects.df[[effect.index,"Effect.Variation"]]
  eff = rnorm(n = 1, mean = real.eff, sd = effvar)
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

# Returns whether or not it is published, based on result and publication bias parameters
published = function(p.value, Alpha, negative.bias) {
  return (p.value <= Alpha | runif(1,0,1) < negative.bias)
}

# Check whether the condition to end the simulation was reached (sample size or number of effects)
reached.sim.end = function(input) {
  if (input$fixed.prev.mode == "none") {
    v = sim.end.tracking(input)
    has_ended = (v >= input$sim.end.value)
  } else {
    ns = sim.end.fixed.tracking(input)
    has_ended = (ns$n_true >= input$sim.end.value * input$fixed.prev &
                   ns$n_false >= input$sim.end.value * (1 - input$fixed.prev))
  }
  has_ended
}

sim.end.fixed.tracking = function (input) {
  if (input$fixed.prev.mode == "above minimum of interest") {
    n_true = estimates.df[Published == T & abs(Real.Effect.Size) >= input$min.effect.of.interest, .N]
    n_false = estimates.df[Published == T & abs(Real.Effect.Size) < input$min.effect.of.interest, .N]
  } else if (input$fixed.prev.mode == "precision") {
    n_true = estimates.df[Published == T & abs(Real.Effect.Size - Estimated.Effect.Size) <= input$min.effect.of.interest, .N]
    n_false = estimates.df[Published == T & abs(Real.Effect.Size - Estimated.Effect.Size) > input$min.effect.of.interest, .N]
  } else if (input$fixed.prev.mode == "non-biased") {
    n_true = estimates.df[Published == T & Biased == F, .N]
    n_false = estimates.df[Published == T & Biased == T, .N]
  } else {
    stop("Invalid value for fixed prevalence mode parameter.")
  }
  list(n_true = n_true, n_false = n_false)
}

sim.end.tracking = function (input) {
  if (input$how.sim.ends == "At a given number of published effects") {
    v = estimates.df[Published == T, .N]
  } else if (input$how.sim.ends == "At a given number of positive published effects") {
    v = estimates.df[Published == T & p.value <= input$alpha.threshold, .N]
  } else if (input$how.sim.ends == "At a given total sample size") {
    v = sum(estimates.df$Sample.Size)
  } else {
    stop("Invalid value for simulation ending criterium.")
  }
  v
}

# Function to simulate the whole experimental procedure, from data generation to analysis and publication
# Calls all the others functions; works as an outline of the main process in the model
perform.experiment = function(effect.index, input) {
  
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
    Sample.Size = input$typical.sample.size,
    Alpha = input$alpha.threshold,
    Min.Interesting.Effect = input$min.effect.of.interest,
    Biased = F,
    Published = published(test$p.value, input$alpha.threshold, input$neg.incentive)
  )
  
  return (new.estimate)
}

# Executes the simulated scientist behavior
scientist.action = function(input) {
  # If simulation ended, return
  if (reached.sim.end(input)) {
    return (1)
  }
  
  # Picks the effect size to be investigated
    effect.index = generate.effect.size(input$sdA, input$weightB, input$meanB, input$sdB, input$interlab.var)

  # Performs the experiment
  xp = perform.experiment(effect.index, input)
  if (input$publish.only.if.large &
      abs(xp$Estimated.Effect.Size) < input$min.effect.of.interest) {
    xp$Published = F
  }
  
  # Bias
  ### With a given probability, if the result is not significant,
  ### change the estimated effect size to be just over the threshold
  if (input$fixed.prev.mode == "non-biased") { target_bias_level = 1 }
  else { target_bias_level = input$bias.level }
  
  if (!xp$Published & runif(1,0,1) < target_bias_level) {
    publishable = F
    while (!publishable) {
      # Redo the experiment until it is significant
      xp = perform.experiment(effect.index, input)
      if (input$publish.only.if.large) {
        publishable = xp$p.value <= input$alpha.threshold & abs(xp$Estimated.Effect.Size) >= input$min.effect.of.interest
      } else {
        publishable = xp$p.value <= input$alpha.threshold
      }
    }
    xp$Biased = T
    xp$Published = T
  }
  
  # Increases the row count, allocates memory if the original allocation is full
  # This is for performance reasons, to avoid growing the data frame (see R Circles of Hell)
  estimates.rowcount <<- estimates.rowcount + 1
  # print(estimates.rowcount)
  if (estimates.rowcount > nrow(estimates.df)) {
    tmp = estimates.df
    tmp[1:nrow(tmp),] = NA
    estimates.df <<- rbindlist(list(estimates.df, tmp), use.names = T, fill = T)
  }

  # Adds new row using data.table::set  
  for (j in 1:ncol(xp)) {
    set(estimates.df, i = as.integer(estimates.rowcount), j, xp[[j]])
  }

  return (0)
}

feedback.message = function(msg, msgtype = "default") {
  print(msg); if (shiny_running) { showNotification(msg, type = msgtype) }
}

# Main function, called from the interface to set the whole model running
run.simulation = function(input) {
  input = sanitize_shiny_input(input)
  
  evdf = data.frame()
  
  input = setup.model(input)
  
  feedback.message("Generating the literature ...")
  
  it = 0
  give.fb = T
  while (it != 1) {
    it = scientist.action(input)
    
    # Progress feedback
    if (input$fixed.prev.mode == "none") {
      progress = round(100 * sim.end.tracking(input) / input$sim.end.value)
    } else {
      ns = sim.end.fixed.tracking(input)
      prog_t = min(c(1, ns$n_true / (input$sim.end.value * input$fixed.prev)))
      prog_f = min(c(1, ns$n_false / (input$sim.end.value * (1 - input$fixed.prev))))
      progress = round(50 * (prog_t + prog_f))
      # print(glue("{prog_t} + {prog_f}"))
    }
    
    if (progress %% 25 == 0 & give.fb) {
      feedback.message(glue("{progress}% done ..."))
      give.fb = F
    } else if (progress %% 25 != 0) {
      give.fb = T
    }
  }

  # Cut non-filled tail of the data table
  estimates.df <<- estimates.df[!is.na(Published),]

  # If the prevalence is fixed, sample accordingly from estimates.df
  if (input$fixed.prev.mode != "none") {
    n_true = round(input$sim.end.value * input$fixed.prev)
    n_false = round(input$sim.end.value * (1 - input$fixed.prev))
    
    if (input$fixed.prev.mode == "above minimum of interest") {
      true.estimates = estimates.df[Published == T & abs(Real.Effect.Size) >= input$min.effect.of.interest,][sample(.N, n_true),]
      false.estimates = estimates.df[Published == T & abs(Real.Effect.Size) < input$min.effect.of.interest,][sample(.N, n_false),]
    } else if (input$fixed.prev.mode == "precision") {
      true.estimates = estimates.df[Published == T & abs(Real.Effect.Size - Estimated.Effect.Size) <= input$repro.detect,][sample(.N, n_true),]
      false.estimates = estimates.df[Published == T & abs(Real.Effect.Size - Estimated.Effect.Size) > input$min.effect.of.interest,][sample(.N, n_false),]
    } else if (input$fixed.prev.mode == "non-biased") {
      true.estimates = estimates.df[Published == T & Biased == F,][sample(.N, n_true),]
      false.estimates = estimates.df[Published == T & Biased == T,][sample(.N, n_false),]
    }
    
    estimates.df <<- rbindlist(list(true.estimates, false.estimates))
  }
  
  # Performs replications
  if (input$calc.repro) {
    feedback.message("Replicating experiments ...")
    replications.df <<- perform.replications(input)
  }
  
  if (shiny_running) {
    feedback.message("Evaluating the literature ...")
    
    new.measures = make.evaluation.tests()
    
    if (input$scenarioName == "") {
      new.measures$scenario = input$scenario 
    } else {
      new.measures$scenario = input$scenarioName
    }
    
    evdf = rbind(evdf, new.measures)
    
    eval.df <<- evdf
  
    if (input$calc.repro) {
      feedback.message("Evaluating replications ...")
      rep.eval.df <<- make.rep.evaluation.tests(input$min.effect.of.interest, input$repro.detect)
    }
  }
  
  feedback.message("... and ... Finished!", "error")
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
