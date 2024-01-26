suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pool))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(ggbeeswarm))


# Override ggplot2 defaults
theme_set(theme_light())
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")

# Connect to the database and establish choices
pool <- dbPool(
  drv = duckdb(),
  dbdir = "appdata/pancan.duckdb", 
  read_only = TRUE
)

onStop(function() { poolClose(pool) })

# Initialize choices
genelist <- readRDS("appdata/genes.rds")
cancers <- readRDS("appdata/cancers.rds")

tbl_choices <- data.frame(
  id = c("expr", "cn", "meth", "mut"),
  names = c("Gene Expression", "Copy Number", "Methylation (1kb from TSS)", "Damaging Mutations")
)

col_choices <- data.frame(
  id = c("N/A", "age", "race", "gender", "ajcc_pathologic_tumor_stage", "vital_status", "cancer_type", "os_time"),
  names = c("N/A", "Age", "Race", "Gender", "Pathologic Stage", "Vital Status", "Cancer Type", "Overall Survival Time")
)


# UI ------------------------------------------------------------------------------------------


ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  titlePanel("TCGA PanCancer Data Explorer (alpha)"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      pickerInput(
        "types",
        label = "Cancer type(s)",
        cancers,
        selected = "BRCA",
        multiple = TRUE,
        options = pickerOptions(actionsBox = TRUE, size = 10)
      ),
      pickerInput(
        "x_gene",
        label = "Gene (x-axis)",
        choices = genelist,
        selected = "CDK12",
        options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE, size = 10)
      ),
      selectInput(
        "x_tbl",
        label = "Data Type (x-axis)",
        choices = setNames(tbl_choices$id, tbl_choices$names),
        selected = "Gene Expression"
      ),
      pickerInput(
        "y_gene",
        label = "Gene (y-axis)",
        choices = genelist,
        selected = "CDK9",
        options = pickerOptions(actionsBox = TRUE, liveSearch = TRUE, size = 10)
      ),
      selectInput(
        "y_tbl",
        label = "Data Type (y-axis)",
        choices = setNames(tbl_choices$id, tbl_choices$names),
        selected = "Gene Expression"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "scatter",
          fluidRow(
            column(
              2,
              numericInput("psize", "Point Size", value = 3, min = 1, max = 10, step = 1),
              numericInput("alpha", "Point Alpha", value = 0.5, min = 0, max = 1, step = 0.1),
              checkboxInput("trend", "Fit Trendline", value = FALSE),
              checkboxInput("smooth", "Fit Smoother", value = FALSE),
              checkboxInput("rug", "Show Rugs", value = FALSE),
              selectInput(
                "scatter_col",
                "Color By",
                choices = setNames(col_choices$id, col_choices$names),
                selected = "N/A"
              ),
              downloadButton("download"),
            ),
            column(10, plotOutput("scatter", brush = "scatter_brush"))
          ),
          gt_output("data")
        ),
        tabPanel(
          "x-axis",
          fluidRow(
            column(
              2,
              numericInput("xsize", "Point Size", value = 1, min = 1, max = 10, step = 1),
              numericInput("xalpha", "Point Alpha", value = 0.5, min = 0.1, max = 1, step = 0.1),
              selectInput(
                "x_col",
                "Color By",
                choices = c("N/A" = "N/A", "Race" = "race", "Gender" = "gender", "Tumor Stage" = "ajcc_pathologic_tumor_stage","Vital Status" = "vital_status"),
                selected = "N/A"
              )
            ),
            column(10, plotOutput("x_axis")),
            gt_output("xdata")
          )
        ),
        tabPanel(
          "y-axis",
          fluidRow(
            column(
              2,
              numericInput("ysize", "Point Size", value = 1, min = 1, max = 10, step = 1),
              numericInput("yalpha", "Point Alpha", value = 0.5, min = 0.1, max = 1, step = 0.1),
              selectInput(
                "y_col",
                "Color By",
                choices = c("N/A" = "N/A", "Race" = "race", "Gender" = "gender", "Tumor Stage" = "ajcc_pathologic_tumor_stage","Vital Status" = "vital_status"),
                selected = "N/A"
              )
            ),
            column(10, plotOutput("y_axis")),
            gt_output("ydata")
          )
        ),
        tabPanel(
          "about",
          htmltools::includeMarkdown("appdata/about.md")
        )
      )
    )
  )
)


# Server --------------------------------------------------------------------------------------


server <- function(input, output, session) {
  data <- reactive({
    df <- inner_join(
      x = {
        tbl(pool, input$x_tbl) |>
          filter(Gene == !!input$x_gene) |>
          semi_join(y = {
            tbl(pool, "clinical") |> filter(cancer_type %in% !!input$types)
          }, by = join_by(SampleBarcode))
      },
      y = {
        tbl(pool, input$y_tbl) |>
          filter(Gene == !!input$y_gene) |>
          semi_join(y = {
            tbl(pool, "clinical") |> filter(cancer_type %in% !!input$types)
          }, by = join_by(SampleBarcode))
      },
      by = join_by(SampleBarcode)
    ) |>
      left_join(y = tbl(pool, "clinical"), by = join_by(SampleBarcode)) |>
      collect()

    validate(
      need(
        nrow(df) != 0,
        "There are no matches in the database. Try selecting a different gene, data type, or adding different samples."
      )
    )

    return(df)
  })

  scatter_plot <- reactive(
    {
      x_var <- switch(input$x_tbl,
        "expr" = "Expression",
        "cn" = "CopyNumber",
        "meth" = "meanBeta",
        "mut" = "Mutations"
      )

      y_var <- switch(input$y_tbl,
        "expr" = "Expression",
        "cn" = "CopyNumber",
        "meth" = "meanBeta",
        "mut" = "Mutations"
      )

      if (input$x_tbl == input$y_tbl) {
        x_var <- paste0(x_var, ".x")
        y_var <- paste0(y_var, ".y")
      }

      p <- ggplot(
        data(),
        aes(
          x = .data[[x_var]],
          y = .data[[y_var]],
          color = if (input$scatter_col == "N/A") NULL else .data[[input$scatter_col]]
        )
      ) +
        geom_point(size = input$psize, pch = 16, alpha = input$alpha) +
        labs(
          title = paste(input$x_gene, "vs.", input$y_gene),
          x = paste0(input$x_gene, " (", tbl_choices$names[tbl_choices$id == input$x_tbl], ")"),
          y = paste0(input$y_gene, " (", tbl_choices$names[tbl_choices$id == input$y_tbl], ")"),
          color = col_choices$names[col_choices$id == input$scatter_col]
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 18),
          plot.subtitle = element_text(size = 12),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "right"
        )

      if (input$trend) {
        mod <- lm(data()[[y_var]] ~ data()[[x_var]])
        b <- round(coef(mod)[1], 2)
        m <- round(coef(mod)[2], 2)
        eq <- paste0("y ~ ", m, "x", " + ", b)
        
        p <- p + geom_smooth(method = "lm", color = "red2") + labs(subtitle = eq)
      }

      if (input$smooth) {
        p <- p + geom_smooth(color = "blue2")
      }
      
      if (input$rug) {
        p <- p + geom_rug(sides = "bl", alpha = 0.1, color = "black")
      }

      return(p)
    }
  )
  output$scatter <- renderPlot(scatter_plot(), res = 100)

  output$data <- render_gt({
    x_var <- switch(input$x_tbl,
      "expr" = "Expression",
      "cn" = "CopyNumber",
      "meth" = "meanBeta",
      "mut" = "Mutations"
    )

    y_var <- switch(input$y_tbl,
      "expr" = "Expression",
      "cn" = "CopyNumber",
      "meth" = "meanBeta",
      "mut" = "Mutations"
    )

    if (input$x_tbl == input$y_tbl) {
      x_var <- paste0(x_var, ".x")
      y_var <- paste0(y_var, ".y")
    }

    keep_cols <- c(
      "SampleBarcode", "Gene.x", x_var, "Gene.y", y_var, "age", "gender", "race",
      "ajcc_pathologic_tumor_stage", "vital_status", "cancer_type", "os_time"
    )
    df <- data()[, keep_cols]

    brushedPoints(
      df,
      input$scatter_brush,
      xvar = x_var,
      yvar = y_var
    ) |>
      gt() |>
      tab_header(
        title = gt::md("**Selected Data**")
      ) |>
      opt_interactive(use_compact_mode = TRUE)
  })

  xplot <- reactive(
    {
      x_var <- switch(input$x_tbl,
        "expr" = "Expression",
        "cn" = "CopyNumber",
        "meth" = "meanBeta",
        "mut" = "Mutations"
      )

      if (input$x_tbl == input$y_tbl) {
        x_var <- paste0(x_var, ".x")
      }

      data() |>
        ggplot(aes(
          x = reorder(cancer_type, .data[[x_var]], median),
          y = .data[[x_var]],
          color = if (input$x_col == "N/A") NULL else .data[[input$x_col]]
        )) +
        geom_quasirandom(size = input$xsize, pch = 16, alpha = input$xalpha) +
        stat_summary(aes(group = cancer_type), fun = median, geom = "crossbar", color = "red2", width = 0.7, lwd = 0.5) +
        labs(
          title = paste("Distribution of", input$x_gene, tbl_choices$names[tbl_choices$id == input$x_tbl]),
          x = "Cancer Type",
          y = tbl_choices$names[tbl_choices$id == input$x_tbl],
          color = col_choices$names[col_choices$id == input$x_col]
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 18),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14)
        )
    }
  )
  output$x_axis <- renderPlot(xplot(), res = 100)
  
  output$xdata <- render_gt({
    x_var <- switch(input$x_tbl,
                    "expr" = "Expression",
                    "cn" = "CopyNumber",
                    "meth" = "meanBeta",
                    "mut" = "Mutations"
    )
    
    if (input$x_tbl == input$y_tbl) {
      x_var <- paste0(x_var, ".x")
    }
    
    if (input$x_col == "N/A" || input$x_col == "cancer_type") {
      df <- data() |> 
        group_by(cancer_type) |> 
        summarise(mean_xvar = mean(.data[[x_var]], na.rm = TRUE),
                  median_xvar = median(.data[[x_var]], na.rm = TRUE),
                  sd_xvar = sd(.data[[x_var]], na.rm = TRUE),
                  min_xvar = min(.data[[x_var]], na.rm = TRUE),
                  max_xvar = max(.data[[x_var]], na.rm = TRUE),
                  N_xvar = n(),
                  .groups = "drop")
    } else {
      df <- data() |> 
        group_by(.data[[input$x_col]], cancer_type) |> 
        summarise(mean_xvar = mean(.data[[x_var]], na.rm = TRUE),
                  median_xvar = median(.data[[x_var]], na.rm = TRUE),
                  sd_xvar = sd(.data[[x_var]], na.rm = TRUE),
                  min_xvar = min(.data[[x_var]], na.rm = TRUE),
                  max_xvar = max(.data[[x_var]], na.rm = TRUE),
                  N_xvar = n(),
                  .groups = "drop")
    }
    
    df |> 
      gt() |> 
      tab_header(
        title = gt::md("**Summary Stats**")
      ) |>
      opt_interactive(use_compact_mode = TRUE)
    
  })

  yplot <- reactive(
    {
      y_var <- switch(input$y_tbl,
        "expr" = "Expression",
        "cn" = "CopyNumber",
        "meth" = "meanBeta",
        "mut" = "Mutations"
      )

      if (input$x_tbl == input$y_tbl) {
        y_var <- paste0(y_var, ".y")
      }

      data() |>
        ggplot(aes(
          x = reorder(cancer_type, .data[[y_var]], median),
          y = .data[[y_var]],
          color = if (input$y_col == "N/A") NULL else .data[[input$y_col]]
        )) +
        geom_quasirandom(size = input$ysize, pch = 16, alpha = input$yalpha) +
        stat_summary(aes(group = cancer_type), fun = median, geom = "crossbar", color = "red2", width = 0.7, lwd = 0.5) +
        labs(
          title = paste("Distribution of", input$y_gene, tbl_choices$names[tbl_choices$id == input$y_tbl]),
          x = "Cancer Type",
          y = tbl_choices$names[tbl_choices$id == input$y_tbl],
          color = col_choices$names[col_choices$id == input$y_col]
        ) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 18),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14)
        )
    }
  )
  output$y_axis <- renderPlot(yplot(), res = 100)
  
  output$ydata <- render_gt({
    y_var <- switch(input$y_tbl,
                    "expr" = "Expression",
                    "cn" = "CopyNumber",
                    "meth" = "meanBeta",
                    "mut" = "Mutations"
    )
    
    if (input$x_tbl == input$y_tbl) {
      y_var <- paste0(y_var, ".x")
    }
    
    if (input$y_col == "N/A" || input$y_col == "cancer_type") {
      df <- data() |> 
        group_by(cancer_type) |> 
        summarise(mean_yvar = mean(.data[[y_var]], na.rm = TRUE),
                  median_yvar = median(.data[[y_var]], na.rm = TRUE),
                  sd_yvar = sd(.data[[y_var]], na.rm = TRUE),
                  min_yvar = min(.data[[y_var]], na.rm = TRUE),
                  max_yvar = max(.data[[y_var]], na.rm = TRUE),
                  N_yvar = n(),
                  .groups = "drop")
    } else {
      df <- data() |> 
        group_by(.data[[input$y_col]], cancer_type) |> 
        summarise(mean_yvar = mean(.data[[y_var]], na.rm = TRUE),
                  median_yvar = median(.data[[y_var]], na.rm = TRUE),
                  sd_yvar = sd(.data[[y_var]], na.rm = TRUE),
                  min_yvar = min(.data[[y_var]], na.rm = TRUE),
                  max_yvar = max(.data[[y_var]], na.rm = TRUE),
                  N_yvar = n(),
                  .groups = "drop")
    }
    
    df |> 
      gt() |> 
      tab_header(
        title = gt::md("**Summary Stats**")
      ) |>
      opt_interactive(use_compact_mode = TRUE)
    
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0("tcga-pancan_", format(Sys.time(), "%Y-%m-%d"), ".zip")
    },
    content = function(file) {
      
      tmp <- tempdir()
      setwd(tmp)
      
      data_file <- data.table::fwrite(data(), "tcga-pancan-data.tsv", sep="\t")
      scatter_plot <- ggsave("scatter-plot.pdf", plot = scatter_plot(), device = "pdf", width = 10, height = 6)
      x_plot <- ggsave("x-axis-violin-plot.pdf", plot = xplot(), device = "pdf", width = 11, height = 7)
      y_plot <- ggsave("y-axis-violin-plot.pdf", plot = yplot(), device = "pdf", width = 11, height = 7)
      files <- c("tcga-pancan-data.tsv", "scatter-plot.pdf", "x-axis-violin-plot.pdf", "y-axis-violin-plot.pdf")
      
      zip(zipfile = file, files = files)
    }, 
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
