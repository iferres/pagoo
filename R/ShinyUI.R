
#' @import shiny
#' @import shinyjs
#' @importFrom DT DTOutput
#' @importFrom plotly plotlyOutput
ui <- fluidPage(
  useShinyjs(),
  # setBackgroundColor(),

  titlePanel("Pagoo Shiny App"),

  tabsetPanel(

    tabPanel(
      "Summary",

      fluidRow(
        column(width = 3, sliderInput("core_level",
                                      "Core Level",
                                      min=85,
                                      max=100,
                                      value=95)),

        column(width = 2,
               radioButtons("input_type",
                            label = "Selection type:",
                            choices = c("organism", "metadata"),
                            selected = "organism")
        ),

        uiOutput("select_organisms")

      ),

      fluidRow(
        column(width = 2,
               align = "center",
               DT::DTOutput("summary_counts"),
               plotlyOutput("core_evol", width = "auto")
        ),
        column(width = 3, offset = 0, plotlyOutput("pie")),
        column(width = 3, offset = 0, plotlyOutput("barplot")),
        column(width = 4, offset = 0, plotlyOutput("curves"))
      ),

      # hr(),

      fluidRow(
        column(width = 6,style = "margin-top:-5em;",
               align = "center",
               DT::DTOutput("core_clusters"), br()
        ),
        column(
          width = 6,style = "margin-top:-5em;",
          align = "center",
          DT::DTOutput("core_genes"), br()
        )
      )


    ), # Ends tabPanel "Summary"

    tabPanel(
      "Accessory",

      fluidRow(
        column(
          width = 3,
          uiOutput(outputId = "accs_slider")
        ),

        column(
          width = 1,
          br(),
          uiOutput("pca_x"),
        ),

        column(
          width = 1,
          br(),
          uiOutput("pca_y")
        ),

        column(
          width = 2,
          # offset = 1,
          # br(),
          checkboxInput("check_color_meta",
                        label = "Color by",
                        value = FALSE),
          uiOutput("pca_color_meta")

        ),

        column(
          width = 4,
          # offset = 1,
          # br(),
          DT::DTOutput("pca_summary")
        )
      ),

      hr(),

      fluidRow(
        splitLayout(
          plotlyOutput("heat_freq", height = "350px"),
          plotlyOutput("pca", height = "350px"),
          cellWidths = c("60%", "40%")
        )
      ),

      hr(),

      fluidRow(
        column(
          width = 6,align = "center",
          DT::DTOutput("accs_clusters"),
        ),
        column(
          width = 6, align = "center",
          DT::DTOutput("accs_genes")
        )

      )

    ) # Ends tabPanel "Analysis"

  ) # Ends tabPanel

)
