# Define UI
ui <- fluidPage(
   
   # Application title
   titlePanel("How Much Reproducibility?"),
   
   theme = "style.css",
   
   tabsetPanel(
     
     tabPanel("Setup & Run",
       
       fluidRow(
         tags$head(
           tags$style(HTML('#runButton{background-color:darkgray; font-weight:bold}'))
         ),
         column(3,
                h4("Underlying Distribution of Effect Sizes"),
                
                plotOutput('EffDistPlot', width = "100%", height = 200),
                # span(textOutput('perc.above.minimum'), style="font-size: 85%"),
                
                sliderInput('sdA', "SD of Distribution A", 
                            min = 0, max = 1, value = 0.1, 
                            step = 0.1, round = -2, ticks = F),
                sliderInput('weightB', "Weight of Distribution B", 
                            min = 0, max = 1, value = 0.5, 
                            step = 0.1, round = -2, ticks = F),
                sliderInput('sdB', "SD of Distribution B", 
                            min = 0, max = 1, value = 0.1, 
                            step = 0.1, round = -2, ticks = F),
                sliderInput('meanB', "Mean of Distribution B", 
                            min = 0, max = 1, value = 1, 
                            step = 0.1, round = -2, ticks = F)
         ),
         
         column(3,
                h4("Power/Sample Size"),
                
                sliderInput('min.effect.of.interest', 'Minimum Effect of Interest',
                            min = 0.1, max = 1, value = 0.5,
                            step = 0.1, round = -2, ticks = F),
                
                sliderInput('typical.power', 'Power',
                            min = 0.1, max = 0.9, value = 0.8, 
                            step = 0.1, round = -1, ticks = F),
                
                span(textOutput('out.N'), style="font-size: 100%"),
                
                br(),
                
                sliderInput('alpha.threshold', 'Alpha', 
                            min = 0.01, max = 0.1, value = 0.05, 
                            step = 0.01, round = -2, ticks = F)
                
         ),
         
         column(3,
                
                h4("Other Parameters"),
                
                sliderInput('interlab.var', 'Interlab Variation',
                            min = 0, max = 1, value = 0, 
                            step = 0.05, round = -2, ticks = F),
                
                sliderInput('bias.level', 'Bias Chance', 
                            min = 0, max = 1, value = 0, 
                            step = 0.05, round = -2, ticks = F),
                
                sliderInput('neg.incentive', 'Negative Results Incentive', 
                            min = 0, max = 1, value = 0, 
                            step = 0.1, round = -1, ticks = F),
                
                hr(),
                
                h4("End Condition"),
                
                selectInput('how.sim.ends', "How should the simulation end?",
                            choices = c("At a given number of published effects",
                                        "At a given total sample size",
                                        "At a given number of positive published effects",
                                        "At a given number of unique positive published effects"),
                            selectize = F),
                
                sliderInput('sim.end.value', 'Set the value to stop the simulation.', 
                            min = 00, max = 5000, value = 200, 
                            step = 100, round = 0, ticks = F)
                
         ),
         
         column(3,
                h4("Setup"),
                
                selectInput('scenario', "Pick a preset scenario (optional)",
                            choices = preset.scenarios, selected = "Small (for quick testing)",
                            selectize = F),
                
                textInput("scenarioName",
                          "Choose a name for your simulation (optional)", ""),
                
                br(),

                checkboxInput('calc.repro', "Calculate reproducibility rates?", value = F),
                sliderInput('repro.repeats', 'Sample size for reproducibility rate:',
                            min = 1, max = 10, value = 1,
                            step = 1, round = 0, ticks = F),
                br(),
                
                actionButton("runButton", "Run simulation!"),
                
                h4("Save & Load"),
                
                downloadButton("downloadData", "Download"),
                
                br(),
                br(),
                
                fileInput("loadDataFile", "Choose ZIP File with results:",
                          multiple = FALSE, placeholder = "", accept = "application/zip")
                
         )
       )
     ),
     
     tabPanel("Results",
        fluidRow(
         
         column(2,
                
                h4("Parameters"),
                
                tableOutput('parameters.table')

         ),
         
         column(3,
                
                h4("Outcomes"),
                
                tableOutput('outcome.table')
         ),
         
         column(5,
                
                h4("Estimated Effect Sizes"),
                
                plotOutput('outcome.plot')
         ),
         
         column(2,
                
                h4("Options"),
                
                checkboxInput('show.published.only', "Show only published results?", value = T),
                checkboxInput('show.biased', "Show biased results?", value = T),
                checkboxInput('show.density', "Show density?", value = F),
                checkboxInput('show.histogram', "Show histogram?", value = T),
                checkboxInput('show.dotplot', "Show dotplot?", value = F),
                br(),
                actionButton("updateButton", "Update Results")
         )
         
       )
     ),
     
     tabPanel("Help",
              sidebarLayout(
                
                sidebarPanel(
                  h4("Overview"),
                  p(
                    HTML(
                      "Aqui teria o texto da ajuda, de forma que as pessoas tenham uma ideia de como o modelo funciona e o que dá pra fazer com ele. <br><br><b>Teste negrito</b> <br> Ao lado elas veriam nas abas a descrição de cada parâmetro e outcome. Parece que o básico de <pre>HTML</pre> funciona bem aqui, então dá pra escrever a ajuda em algum outro editor WYSIWYG e colar aqui.")
                  )
                ),
                
                mainPanel(
                  tabsetPanel(
                    tabPanel("Parameters"),
                    tabPanel("Outcomes")
                  )
                )
              )
            )
        )
)