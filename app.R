library(shiny)
library(plotly)
library(leaflet)
library(ggplot2)
library(MASS)
library(pracma)
library(shinyWidgets)


# Funciones de envoltura
wrap_longitude <- function(lon) {
  ((lon + 180) %% 360) - 180
}

wrap_latitude <- function(lat) {
  ((lat + 90) %% 180) - 90
}

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      ".custom-title {
         font-size: 28px;
         font-weight: bold;
         color: #2c3e50;
         text-align: center;
         margin-bottom: 30px;
         font-family: 'Helvetica Neue', sans-serif;
       }"
    ))
  ),

  div(class = "custom-title", 
      "Trajectory Simulation and Angular Analysis of Multi-Dimensional Integral Fractional Ornstein-Uhlenbeck Process"),

  sidebarLayout(
    sidebarPanel(
      h4("Simulation"),
      numericInput("T", "Total Time (T)", value = 200, min = 1),
      numericInput("n", "Number of Steps (n)", value = 400, min = 10),

      h4("Initial Values"),
      numericInput("mu_ini_1", HTML("&mu;<sub>1</sub>(0)"), value = 0),
      numericInput("mu_ini_2", HTML("&mu;<sub>2</sub>(0)"), value = 0),
      numericInput("mu_ini_3", HTML("&mu;<sub>3</sub>(0)"), value = 0),

      h4("Parameters of Multi-Dimensional Integral Fractional Ornstein-Uhlenbeck Process"),
      numericInput("s1", HTML("&sigma;<sub>1</sub>"), value = 2),
      numericInput("b1", HTML("&beta;<sub>1</sub>"), value = 7),
      numericInput("h1", HTML("H<sub>1</sub>"), value = 0.75),

      numericInput("s2", HTML("&sigma;<sub>2</sub>"), value = 3),
      numericInput("b2", HTML("&beta;<sub>2</sub>"), value = 10),
      numericInput("h2", HTML("H<sub>2</sub>"), value = 0.5),

      numericInput("s3", HTML("&sigma;<sub>3</sub>"), value = 2.5),
      numericInput("b3", HTML("&beta;<sub>3</sub>"), value = 5),
      numericInput("h3", HTML("H<sub>3</sub>"), value = 0.25),

      h4("Cross-Correlation Parameters"),
      numericInput("z12", HTML("z<sub>1,2</sub>"), value = 0.2),
      numericInput("z13", HTML("z<sub>1,3</sub>"), value = -0.5),
      numericInput("z23", HTML("z<sub>2,3</sub>"), value = 0.7),

      actionButton("simulate", "Simulate Trajectory")
    ),

    mainPanel(
      leafletOutput("mapPlot", height = 500),
      plotOutput("anglePlots"),
      uiOutput("correlationTable")
    )
  )
)

# Server logic
server <- function(input, output, session) {

  simulate_data <- eventReactive(input$simulate, {
    source("function.R", local = TRUE)

    T <- input$T
    n <- input$n
    D <- T / n

    telemetry <- simulation_3d(
      mu_ini_1 = input$mu_ini_1, mu_ini_2 = input$mu_ini_2, mu_ini_3 = input$mu_ini_3,
      z_12 = input$z12, z_13 = input$z13, z_23 = input$z23,
      s_1 = input$s1, s_2 = input$s2, s_3 = input$s3,
      b_1 = input$b1, b_2 = input$b2, b_3 = input$b3,
      h_1 = input$h1, h_2 = input$h2, h_3 = input$h3
    )

    mu_1 <- telemetry[1:(n+1)]
    mu_2 <- telemetry[(n+2):(2*n+2)]  
    mu_3 <- telemetry[(2*n+3):(3*n+3)]

    rhos <- coeff(input$z12, input$z13, input$z23, input$h1, input$h2, input$h3)

    list(mu_1 = mu_1, mu_2 = mu_2, mu_3 = mu_3, n = n, rhos = rhos)
  })

  output$mapPlot <- renderLeaflet({
    req(simulate_data())
    data <- simulate_data()

 df <- data.frame(
  lon = wrap_longitude(data$mu_1),
  lat = wrap_latitude(data$mu_2),
  alt = data$mu_3
)

    range_alt <- range(df$alt, na.rm = TRUE)
    if (diff(range_alt) == 0) {
      range_alt <- range_alt + c(-1, 1)
    }
    pal <- colorNumeric(palette = "viridis", domain = range_alt)

    leaflet(df) %>%
      addTiles() %>%
      fitBounds(~min(lon), ~min(lat), ~max(lon), ~max(lat)) %>%
      addPolylines(~lon, ~lat, color = ~pal(alt), weight = 3) %>%
      addCircleMarkers(~lon, ~lat, radius = 3, color = ~pal(alt), fillOpacity = 0.9) %>%
      addLegend("bottomright", pal = pal, values = df$alt, title = "Altitude (m)")
  })

  output$anglePlots <- renderPlot({
    req(simulate_data())
    data <- simulate_data()
    n <- data$n

    mu_1 <- data$mu_1
    mu_2 <- data$mu_2
    mu_3 <- data$mu_3

    angles_xy <- numeric(n)
    angles_xz <- numeric(n)
    angles_yz <- numeric(n)

    for(i in 1:n){
      dx <- mu_1[i+1] - mu_1[i]
      dy <- mu_2[i+1] - mu_2[i]
      dz <- mu_3[i+1] - mu_3[i]

      angles_xy[i] <- atan2(dy, dx)
      angles_xz[i] <- atan2(dz, dx)
      angles_yz[i] <- atan2(dz, dy)
    }

    df_angles <- data.frame(
      angle = c(angles_xy, angles_xz, angles_yz),
      pair = factor(rep(c("XY", "XZ", "YZ"), each = n))
    )

    ggplot(df_angles, aes(x = angle, fill = pair)) +
      geom_histogram(binwidth = pi/18, color = "white") +
      coord_polar(start = 0) +
      facet_wrap(~pair) +
      theme_minimal() +
      labs(title = "Angular Distribution", x = "Angle (radians)", y = "Count")
  })

  output$correlationTable <- renderUI({
    req(simulate_data())
    rhos <- simulate_data()$rhos

    tags$table(style = "width: 60%; margin-top: 20px; border-collapse: collapse;",
               tags$tr(
                 tags$th(colspan = 6, style = "background-color: #d9edf7; padding: 10px; text-align: center; font-size: 16px;",
                         "Cross-Correlation Parameters")
               ),
               tags$tr(
                 tags$th(style = "background-color: #f0f0f0; padding: 8px; text-align: center;", HTML("&rho;<sub>1,2</sub>")),
                 tags$td(style = "padding: 8px; text-align: center;", round(rhos[1], 4)),
                 tags$th(style = "background-color: #f0f0f0; padding: 8px; text-align: center;", HTML("&rho;<sub>1,3</sub>")),
                 tags$td(style = "padding: 8px; text-align: center;", round(rhos[2], 4)),
                 tags$th(style = "background-color: #f0f0f0; padding: 8px; text-align: center;", HTML("&rho;<sub>2,3</sub>")),
                 tags$td(style = "padding: 8px; text-align: center;", round(rhos[3], 4))
               )
    )
  })
}

shinyApp(ui, server)