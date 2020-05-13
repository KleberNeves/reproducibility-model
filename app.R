# Shiny app
source("ui.r", local = TRUE)
source("server.r", local = TRUE)

# Run the application 
shinyApp(ui = ui, server = server, session = session)