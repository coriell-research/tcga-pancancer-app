suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(duckdb))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gt))


con <- dbConnect(duckdb(), "appdata/pancan.duckdb", read_only = TRUE)
genelist <- readRDS("appdata/genes.rds")
cancers <- readRDS("appdata/cancers.rds")


# Define UI for application
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "litera"),
  titlePanel("TCGA PanCancer Data Explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput(
        inputId = "dtype",
        label = "Data Type",
        choices = c(
          "Gene Expression" = "expr",
          "Copy Number" = "cn",
          "Methylation" = "meth"
        ),
        selected = "Gene Expression"
      ),
      selectInput(
        inputId = "cancer",
        label = "Cancer type",
        choices = cancers,
        selected = "BRCA",
        multiple = FALSE
      ),
      selectInput(
        inputId = "xvar",
        label = "Gene (x-axis)",
        choices = genelist,
        selected = "CDK12",
        selectize = TRUE
      ),
      selectInput(
        inputId = "yvar",
        label = "Gene (y-axis)",
        choices = genelist,
        selected = "CDK9",
        selectize = TRUE
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
              selectInput(
                "scatter_col",
                "Color By",
                choices = c(
                  "N/A" = "N/A",
                  "Age" = "age",
                  "Gender" = "gender",
                  "Race" = "race",
                  "Tumor Stage" = "ajcc_pathologic_tumor_stage",
                  "Vital Status" = "vital_status"
                ),
                selected = "N/A"
              )
            ),
            column(10, plotOutput("scatter", brush = "scatter_brush"))
          ),
          gt_output("data")
        ),
        tabPanel(
          "x-axis",
          plotOutput("x_axis")
        ),
        tabPanel(
          "y-axis",
          plotOutput("y_axis")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  data <- reactive({
    colname <- switch(input$dtype,
      expr = "Expression",
      cn = "CopyNumber",
      meth = "meanBeta"
    )

    tbl(con, input$dtype) |>
      filter(Gene %in% c(!!input$xvar, !!input$yvar)) |>
      inner_join(
        y = {
          tbl(con, "clinical") |> filter(cancer_type == !!input$cancer)
        },
        by = "SampleBarcode"
      ) |>
      collect() |>
      pivot_wider(
        id_cols = c(SampleBarcode, age, gender, race, ajcc_pathologic_tumor_stage, vital_status),
        names_from = Gene,
        values_from = sym(colname)
      )
  })

  output$scatter <- renderPlot(
    {
      df <- data()

      p <- ggplot(df, aes(
        x = !!sym(input$xvar),
        y = !!sym(input$yvar),
        color = if (input$scatter_col == "N/A") NULL else !!sym(input$scatter_col)
      )) +
        geom_point(alpha = input$alpha, size = 3) +
        labs(
          title = paste(input$xvar, "vs", input$yvar),
          x = input$xvar,
          y = input$yvar,
          color = input$scatter_color
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 18),
          plot.subtitle = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)
        )

      if (input$trend) {
        fm <- as.formula(paste(input$yvar, "~", input$xvar))
        mod <- lm(fm, data = df)
        b <- round(coef(mod)[1], 2)
        m <- round(coef(mod)[2], 2)

        p <- p +
          geom_smooth(method = "lm", color = "red2") +
          labs(subtitle = paste0("y ~ ", m, "x + ", b))
      }

      if (input$smooth) {
        p <- p + geom_smooth(color = "blue2", method = "loess")
      }

      return(p)
    },
    res = 100
  )

  output$data <- render_gt({
    keep <- brushedPoints(data(), input$scatter_brush)
    selected_samples <- keep[, "SampleBarcode", drop = TRUE]

    tbl(con, "clinical") |>
      filter(SampleBarcode %in% selected_samples) |>
      collect() |>
      gt() |>
      cols_width(
        SampleBarcode ~ px(175),
        patient ~ px(100),
        cancer_type ~ px(150)
      ) |>
      tab_header(
        title = gt::md("**Selected Samples**")
      ) |>
      opt_interactive()
  })

  output$x_axis <- renderPlot({
    colname <- switch(input$dtype,
      expr = "Expression",
      cn = "CopyNumber",
      meth = "meanBeta"
    )

    tbl(con, input$dtype) |>
      filter(Gene == !!input$xvar) |>
      inner_join(y = tbl(con, "clinical"), by = "SampleBarcode") |>
      collect() |>
      ggplot(aes(
        x = reorder(cancer_type, !!sym(colname), median),
        y = !!sym(colname),
        fill = cancer_type
      )) +
      geom_violin(show.legend = FALSE) +
      labs(title = input$xvar, x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, res = 100)
  
  output$y_axis <- renderPlot({
    colname <- switch(input$dtype,
                      expr = "Expression",
                      cn = "CopyNumber",
                      meth = "meanBeta"
    )
    
    tbl(con, input$dtype) |>
      filter(Gene == !!input$yvar) |>
      inner_join(y = tbl(con, "clinical"), by = "SampleBarcode") |>
      collect() |>
      ggplot(aes(
        x = reorder(cancer_type, !!sym(colname), median),
        y = !!sym(colname),
        fill = cancer_type
      )) +
      geom_violin(show.legend = FALSE) +
      labs(title = input$yvar, x = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, res = 100)
  
}

# Run the application
shinyApp(ui = ui, server = server)
