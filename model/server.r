# Define server logic
server <- function(input, output, session) {
  
  # Obtain output data
  outcome.dist = eventReactive(c(input$show.published.only, input$show.biased, input$runButton, input$updateButton), {
    x = estimates.df %>% mutate(`Above Minimum` = abs(Real.Effect.Size) > Min.Interesting.Effect) %>%
      select(1:7, `Above Minimum`, Biased, Published)
    if (input$show.published.only) {
      x = x %>% filter(Published)
    }
    if (!input$show.biased) {
      x = x %>% filter(!Biased)
    }
    x
  })
  
  # Reference underlying distribution
  eff.dist.ref = reactive({
      sample.from.dist(sdA = input$sdA, weightB = input$weightB, meanB = input$meanB,
                       sdB = input$sdB, k = 100000)
  })
 
  # Output plot with estimates
  make.plot = eventReactive(c(input$show.published.only, input$show.density, input$show.biased,
 input$show.dotplot, input$show.histogram, input$runButton, input$updateButton), {
    # browser()
    d = outcome.dist()
    
    if (nrow(d) == 0) { ggplot() }
    else {
      # Significance threshold
      # c_alpha = qt(input$alpha.threshold / 2, df = 2 * input$typical.sample.size - 2) / sqrt(input$typical.sample.size)
      # 
      p = ggplot(d, aes(x = Estimated.Effect.Size))
      
      if (input$show.histogram) {
        p = p + geom_histogram(aes(y = ..density..), bins = 50, color = "black", fill = "lightgray", size = 0.35)
      }
      
      if (input$show.density) {
        p = p + geom_density(size = 0.6)
      }
      
      if (input$show.dotplot) {
        p = p + geom_dotplot(bins = 50, dotsize = 0.6, aes(fill = `Above Minimum`), alpha = 0.75) +
          labs(fill = "Real Above Minimum?")
      }
      
      p = p +
        geom_vline(xintercept = input$min.effect.of.interest, color = "cornflowerblue", linetype = "dashed", size = 1, alpha = 0.9) +
        geom_vline(xintercept = -input$min.effect.of.interest, color = "cornflowerblue", linetype = "dashed", size = 1, alpha = 0.9) +
        # geom_vline(xintercept = c_alpha, color = "firebrick", linetype = "dashed", size = 1, alpha = 0.9) +
        # geom_vline(xintercept = -c_alpha, color = "firebrick", linetype = "dashed", size = 1, alpha = 0.9) +
        labs(x = "Effect Size", y = "Frequency") +
        theme_linedraw()
      
      # if (param.df[param.df$Parameter == "effects.dist","Value"] %in% c("1","2","3","5")) {
        p = p + geom_density(data = data.frame(x = eff.dist.ref()), aes(x = x), size = 1)
      # }
      p
    }
  })
  
  # Input underlying distribution plot 
  make.underlying.plot = eventReactive(c(input$sdA, input$sdB, input$meanB, input$weightB, input$min.effect.of.interest), {
    df = data.frame(x = eff.dist.ref())
    am = mean(abs(df$x) >= input$min.effect.of.interest, na.rm = T)
    p = ggplot(df, aes(x = x)) + geom_density() + theme_linedraw() +
      geom_vline(xintercept = -input$min.effect.of.interest,
                 linetype = "dashed", color = "blue") +
      geom_vline(xintercept = input$min.effect.of.interest,
                 linetype = "dashed", color = "blue") +
      labs(x = "Effect Size", y = "",
           title = paste0("Effects above minimum: ", as.character(round(100 * am)), "%"))
    p
  })
  
  # Outcomes table
  make.table = eventReactive(c(input$show.published.only, input$runButton, input$loadDataFile, input$updateButton), {
    browser()
    
    target = c("", "PPVMean", "NPVMean")#, "TypeIErrorRateMean", "TypeIIErrorRateMean")
    
    x = eval.df[eval.df$variable %in% target, ]
    x[x$variable == "PPVMean", "variable"] = "PPV"
    x[x$variable == "NPVMean", "variable"] = "NPV"
    x[x$variable == "TypeIErrorRateMean", "variable"] = "TypeIErrorRate"
    x[x$variable == "TypeIIErrorRateMean", "variable"] = "TypeIIErrorRate"
    x$Measure = gsub("Reproducibility Rate ", "", x$Measure)
    x$Measure[x$variable != ""] = paste(x$variable[x$variable != ""], x$Measure[x$variable != ""])
    
    x = x %>% dcast(Measure ~ Published.Effects.Only., value.var = "Value")
    x = x %>% select(1,3,2) %>% rename(Published = `TRUE`, All = `FALSE`)

    target = c("True Positives", "True Negatives", "False Positives", "False Negatives", "Signal Errors", "Exaggeration Factor", "Signal Errors (Effects > Min)", "Exaggeration Factor (Effects > Min)", "Positive Predictive Value", "Negative Predictive Value", "Discovery Effectiveness")#,"PPV (BRI)", "PPV (BRI) All Papers", "PPV (SSS)", "PPV (SSS) All Papers", "PPV (ST)", "PPV (ST) All Papers", "NPV (BRI)", "NPV (BRI) All Papers", "NPV (SSS)", "NPV (SSS) All Papers", "NPV (ST)", "NPV (ST) All Papers", "TypeIErrorRate (BRI)", "TypeIErrorRate (BRI) All Papers", "TypeIErrorRate (SSS)", "TypeIErrorRate (SSS) All Papers", "TypeIErrorRate (ST)", "TypeIErrorRate (ST) All Papers", "TypeIIErrorRate (BRI)", "TypeIIErrorRate (BRI) All Papers", "TypeIIErrorRate (SSS)", "TypeIIErrorRate (SSS) All Papers", "TypeIIErrorRate (ST)", "TypeIIErrorRate (ST) All Papers") 

    x[2:3] = round(x[2:3], 2)
    mt = match(target, x$Measure)
    x = x[mt[!is.na(mt)],]
    
    if (input$calc.repro) {
      x2 = eval.df[eval.df$variable != "",]
      x2 = x2 %>% dcast(Measure ~ variable, value.var = "Value")
      x2$Text = paste(round(x2$MeanValue, 2), " [", round(x2$CI.low, 2), ", ", round(x2$CI.high, 2), "]", sep = "")
      x2 = x2 %>% select(Measure, Text) %>% rename(Published = Text) %>% mutate(All = NA)
      
      x = rbind(x,x2)
    }
    
    x
  })
  
  # Parameters table
  make.param.table = eventReactive(c(input$runButton, input$updateButton), {
    x = param.df
    # x = x %>% select(2,1)
    x$Parameter = as.character(x$Parameter)
    x$Parameter = as.character(car::recode(x$Parameter, "'n.scientists' = '---'; 'neg.incentive' = 'Negative Publication Incentive'; 'runButton' = '---'; 'measure.error' = 'Measurement Error'; 'rep.incentive' = 'Replication Incentive'; 'loadFilename' = '---'; 'scenarioName' = '---'; 'min.effect.of.interest' = 'Minimum Effect of Interest'; 'sim.end.value' = 'Simulation Stop Value'; 'show.histogram' = '---'; 'dist.param.1' = 'dist.param.1'; 'dist.param.2' = 'dist.param.2'; 'scenario' = '---'; 'show.published.only' = '---'; 'typical.sample.size' = 'Sample Size'; 'typical.power' = 'Power'; 'interlab.var' = 'Interlab Variation'; 'downloadData' = '---'; 'show.density' = '---'; 'alpha.threshold' = 'Alpha'; 'how.sim.ends' = 'How the Simulation Ends'; 'show.dotplot' = '---'; 'bias.level' = 'Bias'; 'uploadData' = '---'; 'saveFilename' = '---'; 'sdA' = 'SD of Distribution A'; 'sdB' = 'SD of Distribution B'; 'meanB' = 'Mean of Distribution B'; 'weightB' = 'Weight of Distribution B'"))
    # browser()
    target = c("SD of Distribution A", "Weight of Distribution B", "SD of Distribution B",
               "Mean of Distribution B", "Minimum Effect of Interest", "% Above minimum",
               "Power", "Sample Size", "Alpha",
               "Negative Publication Incentive", "Replication Incentive",
               "Bias", "Interlab Variation", "Simulation Stop Value", "How the Simulation Ends"
               )
    x = x[which(x$Parameter %in% target),]
    x = x[match(target, x$Parameter),]

    x
  })
  
  # Updating N based on power
  output$out.N = renderText({
    
    # TODO Similar code exists in setup.model, should abstract this into a function
    dist.n = 10000
    a.sample = abs(sample.from.dist(
      input$sdA, input$weightB, input$meanB, input$sdB, dist.n))
    
    # Mean effect size for the effects above the minimum of interest
    mes = mean(a.sample[a.sample >= input$min.effect.of.interest])
    if (is.nan(mes)) { mes = input$min.effect.of.interest }
    
    new = ceiling(power.t.test(
      n = NULL, delta = mes,
      sd = 1, sig.level = input$alpha.threshold,
      power = input$typical.power
    )$n)
    
    if (is.numeric(new)) {
      paste("Sample size (per group): ", new, sep = "")
    } else {
      paste("Sample size (per group): ???", sep = "")
    }
  })
  
  # Save simulation results
  output$downloadData = downloadHandler(
    filename = "results.zip",
    content = function(fname) {
      fs = c()
      tmpdir = tempdir()
      fn = paste(tmpdir,"/estimates.csv", sep = ""); fs = c(fs, fn)
      write.table(x = presynthesis.df, file = fn, sep = ";", row.names = F)
      fn = paste(tmpdir,"/eval.csv", sep = ""); fs = c(fs, fn)
      write.table(x = eval.df, file = fn, sep = ";", row.names = F)
      fn = paste(tmpdir,"/pars.csv", sep = ""); fs = c(fs, fn)
      write.table(x = param.df, file = fn, sep = ";", row.names = F)
      
      t = "Saved!"
      print(t); if (shiny_running) { showNotification(t, type = "default") }
      
      zipr(zipfile = fname, files = fs)
    },
    contentType = "application/zip"
  )
  
  # Load simulation results and update output correspondingly
  load.results = function() {
    inFile = input$loadDataFile
    
    if (is.null(inFile))
      return(NULL)
    
    tmpdir = tempdir()
    unzip(zipfile = inFile$datapath, exdir = tmpdir)

    presynthesis.df <<- read.table(
      file = paste(tmpdir, "/estimates.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
    eval.df <<- read.table(
      file = paste(tmpdir, "/eval.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
    param.df <<- read.table(
      file = paste(tmpdir, "/pars.csv", sep = ""), sep = ";", stringsAsFactors = F, header = T)
    
    t = "Synthesizing ..."
    print(t); if (shiny_running) { showNotification(t, type = "default") }
    
    estimates.df <<- synthesize(presynthesis.df)
    
    t = "Loaded!"
    print(t); if (shiny_running) { showNotification(t, type = "default") }
  }
  
  # Update and interact with outputs
  output$typical.power = renderText({
    pwr = power.t.test(
      n = input$typical.sample.size, delta = input$min.effect.of.interest,
      sd = 1, sig.level = input$alpha.threshold, power = NULL
    )
    if (is.numeric(pwr$power)) {
      paste("Power (sampling error only): ", round(pwr$power * 100, 0), "%", sep = "")
    } else {
      paste("Power (sampling error only): ???", sep = "")
    }
  })
  
  output$typical.power2 = renderText({
    pwr = power.t.test(
      n = input$typical.sample.size, delta = input$min.effect.of.interest,
      sd = (1 + input$measure.error^2 + input$interlab.var^2)^0.5, sig.level = input$alpha.threshold, power = NULL
    )
    if (is.numeric(pwr$power)) {
      paste("Power (all error types): ", round(pwr$power * 100, 0), "%", sep = "")
    } else {
      paste("Power (all error types): ???", sep = "")
    }
  })
  
  observeEvent(input$runButton, {
    run.simulation(input)
  })
  
  observeEvent(input$loadDataFile, {
    load.results()
  })
  
  output$outcome.plot = renderPlot({
    make.plot()
  })
  
  output$EffDistPlot = renderPlot({
    make.underlying.plot()
  })
  
  output$outcome.table = renderTable({
    make.table()
  })
  
  output$parameters.table = renderTable({
    make.param.table()
  })

  
}