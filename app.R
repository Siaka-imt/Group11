# app.R - Final Project Dashboard
# Comparison of Normal vs. Randomized Matrix Approximation

library(shiny)
library(shinydashboard)
library(shinythemes)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Load helper modules
source("algorithms.R")
source("data_utils.R")
source("viz_utils.R")

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Matrix Approx Lab", titleWidth = 300),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Performance Comparison", tabName = "compare", icon = icon("balance-scale")),
      menuItem("Deep Dive (Analysis)", tabName = "deep_dive", icon = icon("microscope"))
    ),
    
    div(style = "padding: 15px;",
        h4("Control Panel"),
        
        # 1. Dataset Selection (Updated for Real Data)
        selectInput("dataset", "1. Select Dataset:",
                    choices = c("MovieLens 100k (Real)", "Volcano Topography (Real)")),
        
        # 2. Algorithm Family Selection (Updated terminology to 'Normal')
        selectInput("method_family", "2. Algorithm Family:",
                    choices = c("SVD (Normal vs Randomized)" = "svd", 
                                "NMF (Random vs rSVD Init)" = "nmf",
                                "CUR (Top vs Weighted)" = "cur")),
        
        # 3. Rank Slider
        sliderInput("rank", "3. Target Rank (k):", min = 2, max = 30, value = 10, step = 1),
        
        hr(),
        p(class = "text-muted", "Click below to run both Normal and Randomized versions."),
        
        actionButton("run", "Run Benchmark", icon = icon("play"), 
                     class = "btn btn-primary btn-block", style = "margin-top: 10px;"),
        
        hr(),
        strong("Status:"),
        verbatimTextOutput("status_info")
    )
  ),
  
  dashboardBody(
    theme = shinytheme("flatly"),
    
    tabItems(
      # ==== TAB 1: ABOUT ====
      tabItem(tabName = "about",
              fluidRow(
                box(width = 12, status = "primary", solidHeader = TRUE,
                    title = "Project Overview",
                    h3("Low-Rank Matrix Approximation"),
                    p("This dashboard implements and compares deterministic vs. randomized algorithms for matrix decomposition."),
                    h4("Implemented Algorithms"),
                    tags$ul(
                      tags$li(strong("SVD:"), "Standard Normal vs. Randomized SVD"),
                      tags$li(strong("NMF:"), "Standard Random Init vs. Optimized (rSVD) Init"),
                      tags$li(strong("CUR:"), "Deterministic (Top Leverage) vs. Randomized Sampling")
                    ),
                    h4("Data Source"),
                    verbatimTextOutput("data_info")
                )
              )
      ),
      
      # ==== TAB 2: COMPARISON (HEAD-TO-HEAD) ====
      tabItem(tabName = "compare",
              fluidRow(
                valueBoxOutput("cmp_time_winner", width = 6),
                valueBoxOutput("cmp_error_diff", width = 6)
              ),
              fluidRow(
                # Normal Column
                box(width = 6, status = "info", solidHeader = TRUE, 
                    title = "Normal / Deterministic",
                    h4(textOutput("lbl_norm_name")),
                    tableOutput("tbl_norm_metrics"),
                    plotOutput("plot_cmp_norm", height = "350px")
                ),
                # Randomized Column
                box(width = 6, status = "warning", solidHeader = TRUE, 
                    title = "Randomized / Optimized",
                    h4(textOutput("lbl_rand_name")),
                    tableOutput("tbl_rand_metrics"),
                    plotOutput("plot_cmp_rand", height = "350px")
                )
              )
      ),
      
      # ==== TAB 3: DEEP DIVE (SHOWING NORMAL / DETERMINISTIC RESULT) ====
      tabItem(tabName = "deep_dive",
              h4("Detailed Analysis of Normal (Deterministic) Algorithm"),
              fluidRow(
                valueBoxOutput("dd_time", width = 4),
                valueBoxOutput("dd_error", width = 4),
                valueBoxOutput("dd_dims", width = 4)
              ),
              fluidRow(
                box(width = 6, status = "info", solidHeader = TRUE,
                    title = "Original Matrix",
                    plotOutput("plot_dd_orig", height = "400px")
                ),
                box(width = 6, status = "info", solidHeader = TRUE,
                    title = textOutput("txt_dd_recon_title"),
                    plotOutput("plot_dd_recon", height = "400px")
                )
              ),
              fluidRow(
                box(width = 6, status = "success", solidHeader = TRUE,
                    title = "Factor Components / Basis (Deterministic)",
                    plotOutput("plot_dd_comp", height = "400px")
                ),
                box(width = 6, status = "warning", solidHeader = TRUE,
                    title = "Singular Values / Residuals",
                    plotOutput("plot_dd_behavior", height = "400px")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive storage for results
  store <- reactiveValues(
    data_list = NULL,
    res_norm = NULL,
    res_rand = NULL,
    time_norm = 0,
    time_rand = 0,
    has_run = FALSE
  )
  
  output$status_info <- renderText({
    if(store$has_run) "Done." else "Ready."
  })
  
  # ==== MAIN EXECUTION BLOCK ====
  observeEvent(input$run, {
    req(input$dataset, input$method_family)
    
    withProgress(message = 'Running Algorithms...', value = 0, {
      
      # 1. Load Data
      incProgress(0.1, detail = "Loading dataset...")
      store$data_list <- load_data(input$dataset)
      X <- store$data_list$X
      k <- input$rank
      
      # 2. Run Normal Version (Deterministic)
      incProgress(0.3, detail = "Running Normal version...")
      t1 <- Sys.time()
      if (input$method_family == "svd") {
        store$res_norm <- algo_svd_normal(X, k)
      } else if (input$method_family == "nmf") {
        store$res_norm <- algo_nmf_normal(X, k)
      } else if (input$method_family == "cur") {
        store$res_norm <- algo_cur_normal(X, k)
      }
      store$time_norm <- as.numeric(Sys.time() - t1)
      
      # 3. Run Randomized Version
      incProgress(0.6, detail = "Running Randomized version...")
      t2 <- Sys.time()
      if (input$method_family == "svd") {
        store$res_rand <- algo_svd_randomized(X, k)
      } else if (input$method_family == "nmf") {
        store$res_rand <- algo_nmf_randomized(X, k)
      } else if (input$method_family == "cur") {
        store$res_rand <- algo_cur_randomized(X, k)
      }
      store$time_rand <- as.numeric(Sys.time() - t2)
      
      store$has_run <- TRUE
      incProgress(1.0, detail = "Done!")
      
      showNotification("Decomposition Complete!", type = "message")
    })
  })
  
  # ==== TAB 1: ABOUT ====
  output$data_info <- renderPrint({
    req(store$data_list)
    cat(store$data_list$description)
  })
  
  # ==== TAB 2: COMPARISON LOGIC ====
  
  output$lbl_norm_name <- renderText({ req(store$res_norm); store$res_norm$name })
  output$lbl_rand_name <- renderText({ req(store$res_rand); store$res_rand$name })
  
  # Metrics Table Helper
  make_metrics <- function(res, time, X) {
    err <- norm(X - res$recon, type = "F") / norm(X, type = "F")
    data.frame(
      Metric = c("Time (s)", "Rel. Error"),
      Value = c(sprintf("%.4f", time), sprintf("%.4f", err))
    )
  }
  
  output$tbl_norm_metrics <- renderTable({
    req(store$res_norm)
    make_metrics(store$res_norm, store$time_norm, store$data_list$X)
  }, colnames = FALSE)
  
  output$tbl_rand_metrics <- renderTable({
    req(store$res_rand)
    make_metrics(store$res_rand, store$time_rand, store$data_list$X)
  }, colnames = FALSE)
  
  output$plot_cmp_norm <- renderPlot({
    req(store$res_norm)
    plot_matrix_heatmap(store$res_norm$recon, "Normal Reconstruction")
  })
  
  output$plot_cmp_rand <- renderPlot({
    req(store$res_rand)
    plot_matrix_heatmap(store$res_rand$recon, "Randomized Reconstruction")
  })
  
  # Winner Value Boxes
  output$cmp_time_winner <- renderValueBox({
    req(store$has_run)
    ratio <- store$time_norm / store$time_rand
    txt <- if(ratio > 1) paste0("Randomized is ", round(ratio, 1), "x Faster") else "Normal is Faster"
    color <- if(ratio > 1) "green" else "yellow"
    valueBox(txt, "Speed Comparison", icon = icon("tachometer-alt"), color = color)
  })
  
  output$cmp_error_diff <- renderValueBox({
    req(store$has_run)
    err_n <- norm(store$data_list$X - store$res_norm$recon, "F")
    err_r <- norm(store$data_list$X - store$res_rand$recon, "F")
    
    # Check if errors are roughly equal (within 5%)
    if (abs(err_n - err_r)/err_n < 0.05) {
      txt <- "Similar Accuracy"
      color <- "blue"
    } else if (err_r < err_n) {
      txt <- "Randomized is More Accurate"
      color <- "green"
    } else {
      txt <- "Normal is More Accurate"
      color <- "yellow"
    }
    valueBox(txt, "Accuracy Comparison", icon = icon("check-circle"), color = color)
  })
  
  # ==== TAB 3: DEEP DIVE LOGIC (SWITCHED TO NORMAL / DETERMINISTIC) ====
  
  output$dd_time <- renderValueBox({
    req(store$time_norm) # Using Normal Time
    valueBox(sprintf("%.3f s", store$time_norm), "Compute Time", icon("clock"), color = "aqua")
  })
  
  output$dd_error <- renderValueBox({
    req(store$res_norm) # Using Normal Result
    err <- norm(store$data_list$X - store$res_norm$recon, "F") / norm(store$data_list$X, "F")
    valueBox(sprintf("%.4f", err), "Relative Error", icon("chart-line"), color = "purple")
  })
  
  output$dd_dims <- renderValueBox({
    req(store$res_norm) # Using Normal Result
    dim_str <- paste(dim(store$res_norm$components), collapse = " x ")
    valueBox(dim_str, "Feature Dims", icon("ruler-combined"), color = "teal")
  })
  
  output$plot_dd_orig <- renderPlot({
    req(store$data_list)
    plot_matrix_heatmap(store$data_list$X, "Original Data Matrix")
  })
  
  output$plot_dd_recon <- renderPlot({
    req(store$res_norm) # Plotting Normal Reconstruction
    plot_matrix_heatmap(store$res_norm$recon, "Reconstructed (Normal/Deterministic)")
  })
  
  output$txt_dd_recon_title <- renderText({
    paste("Reconstructed Matrix (Rank", input$rank, ") - Normal")
  })
  
  output$plot_dd_comp <- renderPlot({
    req(store$res_norm) # Plotting Normal Components
    comps <- store$res_norm$components
    if (ncol(comps) > 10) comps <- comps[, 1:10, drop=FALSE]
    plot_matrix_heatmap(comps, "Basis Components (Top 10)")
  })
  
  output$plot_dd_behavior <- renderPlot({
    req(store$res_norm) # Using Normal Data for Behavior Analysis
    
    if (input$method_family == "svd") {
      # For Normal SVD, plotting the true Singular Values (Scree Plot)
      # We recompute exact singular values for visualization purposes
      s_vals <- svd(store$data_list$X, nu=0, nv=0)$d[1:min(50, length(store$data_list$X))]
      df <- data.frame(Index = 1:length(s_vals), Value = s_vals)
      ggplot(df, aes(Index, Value)) + geom_col(fill="blue") + 
        labs(title = "True Singular Values (Scree Plot)", y="Value") + theme_minimal()
      
    } else if (input$method_family == "nmf") {
      # For Normal NMF, plot Residual Distribution
      residuals <- as.vector(store$data_list$X - store$res_norm$recon)
      if(length(residuals) > 5000) residuals <- sample(residuals, 5000)
      df <- data.frame(Resid = residuals)
      ggplot(df, aes(x=Resid)) + geom_histogram(fill="steelblue", bins=30) +
        labs(title = "Residual Distribution (Normal NMF)", x="Error Value") + theme_minimal()
      
    } else {
      # For Normal CUR, plot Leverage Scores
      levs <- rowSums(store$res_norm$components^2) 
      df <- data.frame(Index = 1:length(levs), Leverage = levs)
      ggplot(df, aes(Index, Leverage)) + geom_point(alpha=0.5) +
        labs(title = "Row Leverage Scores (Deterministic Selection)", y="Score") + theme_minimal()
    }
  })
}

shinyApp(ui, server)