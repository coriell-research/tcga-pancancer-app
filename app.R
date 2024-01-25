suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gt))
suppressPackageStartupMessages(library(ggbeeswarm))


# Override ggplot2 defaults
theme_set(theme_minimal())
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")

scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = as.character(paletteer::paletteer_d("ggsci::default_igv")))
}
scale_color_discrete <- function(...) {
  scale_color_manual(..., values = as.character(paletteer::paletteer_d("ggsci::default_igv")))
}

# Connect to the database and establish choices
con <- dbConnect(duckdb(), "appdata/pancan.duckdb", read_only = TRUE)
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
  titlePanel("TCGA PanCancer Data Explorer"),
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
          "Scatter",
          fluidRow(
            column(
              2,
              numericInput("alpha", "Point Alpha", value = 0.5, min = 0, max = 1, step = 0.1),
              checkboxInput("trend", "Fit Trendline", value = FALSE),
              checkboxInput("smooth", "Fit Smoother", value = FALSE),
              checkboxInput("rug", "Show Rugs", value = FALSE),
              selectInput(
                "scatter_col",
                "Color By",
                choices = setNames(col_choices$id, col_choices$names),
                selected = "N/A"
              )
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
        tbl(con, input$x_tbl) |>
          filter(Gene == !!input$x_gene) |>
          semi_join(y = {
            tbl(con, "clinical") |> filter(cancer_type %in% !!input$types)
          }, by = join_by(SampleBarcode))
      },
      y = {
        tbl(con, input$y_tbl) |>
          filter(Gene == !!input$y_gene) |>
          semi_join(y = {
            tbl(con, "clinical") |> filter(cancer_type %in% !!input$types)
          }, by = join_by(SampleBarcode))
      },
      by = join_by(SampleBarcode)
    ) |>
      left_join(y = tbl(con, "clinical"), by = join_by(SampleBarcode)) |>
      collect()

    validate(
      need(
        nrow(df) != 0,
        "There are no matches in the database. Try selecting a different gene, data type, or adding different samples."
      )
    )

    return(df)
  })


  output$scatter <- renderPlot(
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
        geom_point(size = 3, pch = 16, alpha = input$alpha) +
        labs(
          title = paste(input$x_gene, "vs.", input$y_gene),
          x = paste0(input$x_gene, " (", tbl_choices$names[tbl_choices$id == input$x_tbl], ")"),
          y = paste0(input$y_gene, " (", tbl_choices$names[tbl_choices$id == input$y_tbl], ")"),
          color = col_choices$names[col_choices$id == input$scatter_col]
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "right"
        )

      if (input$trend) {
        p <- p + geom_smooth(method = "lm", color = "red2")
      }

      if (input$smooth) {
        p <- p + geom_smooth(method = "gam")
      }
      
      if (input$rug) {
        p <- p + geom_rug(sides = "bl", alpha = 0.1, color = "black")
      }

      return(p)
    },
    res = 100
  )


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


  output$x_axis <- renderPlot(
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
        geom_quasirandom(size = 3, pch = 16, alpha = 0.8) +
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
    },
    res = 100
  )
  
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
  

  output$y_axis <- renderPlot(
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
        geom_quasirandom(size = 3, pch = 16, alpha = 0.8) +
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
    },
    res = 100
  )
  
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
}

# Run the application
shinyApp(ui = ui, server = server)
