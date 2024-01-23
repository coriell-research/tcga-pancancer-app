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
        inputId = "types",
        label = "Cancer type(s)",
        choices = cancers,
        selected = "BRCA",
        multiple = TRUE
      ),
      selectInput(
        inputId = "x_tbl",
        label = "Data Type (x-axis)",
        choices = c(
          "Gene Expression" = "expr",
          "Copy Number" = "cn",
          "Methylation" = "meth"
        ),
        selected = "Gene Expression"
      ),
      selectizeInput(
        inputId = "xvar",
        label = "Gene (x-axis)",
        choices = genelist,
        selected = "CDK12",
        options = NULL
      ),
      selectInput(
        inputId = "y_tbl",
        label = "Data Type (y-axis)",
        choices = c(
          "Gene Expression" = "expr",
          "Copy Number" = "cn",
          "Methylation" = "meth"
        ),
        selected = "Gene Expression"
      ),
      selectizeInput(
        inputId = "yvar",
        label = "Gene (y-axis)",
        choices = genelist,
        selected = "CDK9"
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
  # data <- reactive()
  # output$scatter <- renderPlot()
  # output$data <- render_gt()
  # output$x_axis <- renderPlot()
  # output$y_axis <- renderPlot()
}

# Run the application
shinyApp(ui = ui, server = server)
