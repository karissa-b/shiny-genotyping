library(shiny)
library(tidyverse)
library(vroom)
library(scales)
library(shinythemes)

# user interface ----------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(

  # theme of the app
  theme = bslib::bs_theme(bootswatch = "simplex"),

  # Application title
  titlePanel("Melt curve genotyping (naglu)",
             windowTitle = "Genotyping from melt curves"
             ),


# tabset input ------------------------------------------------------------

  tabsetPanel(
    tabPanel(title = "Input data",
             fileInput("files",
                       label = "Choose TXT files",
                       multiple = TRUE,
                       accept = c("text/plain", ".txt")),
           tableOutput("files")
    ),

# tabset Inputparams ------------------------------------------------------------
    tabPanel(title = "Params",

             sidebarLayout(
               sidebarPanel(
             # mut temp min
             numericInput(inputId = "mut_temp_min",
                          label = "Mutant peak min temp",
                          value = 80.2,
                          min = 60,
                          max = 100),

             # mut temp max
             numericInput(inputId = "mut_temp_max",
                          label = "Mutant peak max temp",
                          value = 81.9,
                          min = 60,
                          max = 100),

             # wt temp min
             numericInput(inputId = "wt_temp_min",
                          label = "wt peak min temp",
                          value = 82.8,
                          min = 60,
                          max = 100),

             # wt temp max
             numericInput(inputId = "wt_temp_max",
                          label = "wt peak max temp",
                          value = 84.5,
                          min = 60,
                          max = 100),

             # mut deriv min
             numericInput(inputId = "mut_deriv_min",
                          label = "Mutant peak min deriv",
                          value = 27000,
                          min = 0,
                          max = 500000),

             # mut deriv max
             numericInput(inputId = "mut_deriv_max",
                          label = "Mutant peak max deriv",
                          value = 110000,
                          min = 0,
                          max = 500000),

             # wt deriv min
             numericInput(inputId = "wt_deriv_min",
                          label = "wt peak min deriv",
                          value = 40000,
                          min = 0,
                          max = 500000),

             # wt deriv max
             numericInput(inputId = "wt_deriv_max",
                          label = "wt peak max deriv",
                          value = 135000,
                          min = 0,
                          max = 500000)
               ),

             mainPanel(plotOutput("myPlot", width = "100%")))),


# tabset call genotypes ---------------------------------------------------

tabPanel(
  title = "Called genotypes",
  plotOutput("genotypeplot", width = "100%"),
  downloadButton(
    outputId = "downloadGeno",
    label = "Download Genotype table")

)
    ))


    # server function ---------------------------------------------------------


# Define server logic required to draw a histogram
server <- function(input, output, session) {

  output$files <- renderTable(input$files)

  # Reactive expression to read and process the selected files
  data <- reactive({
    req(input$files)  # Ensure files are selected

    # Get the file paths 9which are weird temp files now) and OG names
    files <- input$files

    # check whether I am importing the right thing
    # print(files)

    # Read and process the files using your existing code
    data <- read_delim(
      files %>%
        dplyr::filter(grepl(name, pattern = "deriv")) %>%
        .$datapath,
      delim = "\t",
      skip = 8
    ) %>%
      pivot_longer(
        names_to = "Reading",
        values_to = "deriv",
        starts_with("Reading")
      ) %>%
      dplyr::select(
        position = `Well Location`,
        Reading,
        deriv
      ) %>%
      left_join(
        read_delim(
          files %>%
            dplyr::filter(grepl(name, pattern = "temper")) %>%
            .$datapath,
          delim = "\t",
          skip = 8
        ) %>%
          pivot_longer(
            names_to = "Reading",
            values_to = "temp",
            starts_with("Reading")
          ) %>%
          dplyr::select(
            position = `Well Location`,
            Reading,
            temp
          )
      ) %>%
      left_join(
        read_delim(
          files %>%
            dplyr::filter(grepl(name, pattern = "normal")) %>%
            .$datapath,
          delim = "\t",
          skip = 8
        ) %>%
          pivot_longer(
            names_to = "Reading",
            values_to = "fluor",
            starts_with("Reading")
          ) %>%
          dplyr::select(
            position = `Well Location`,
            Reading,
            fluor
          )
      )
  })

  output$myPlot <- renderPlot({

    # Access reactive data
    data <- data()

    # access reactive params
    mut_temp_min <- input$mut_temp_min
    mut_temp_max <- input$mut_temp_max

    wt_temp_max <- input$wt_temp_max
    wt_temp_min <- input$wt_temp_min

    wt_deriv_max <- input$wt_deriv_max
    wt_deriv_min <- input$wt_deriv_min

    mut_deriv_max <- input$mut_deriv_max
    mut_deriv_min <- input$mut_deriv_min

    # Generate the plot
    ggp <- data %>%
      ggplot(aes(x = temp, y = deriv)) +
      geom_line(aes(group = position)) +
      scale_y_continuous(labels = comma,
                         breaks = seq(0, 500000, by = 25000)) +
      scale_x_continuous(limits = c(78, 87),
                         breaks = seq(1:90)) +
      # mut allele
      annotate(
        geom = "rect",
        alpha = 0.5,
        fill = "red",
        xmin = mut_temp_min,
        xmax = mut_temp_max,
        ymin = mut_deriv_min,
        ymax = mut_deriv_max
      ) +
      # wt allele
      annotate(
        geom = "rect",
        alpha = 0.5,
        fill = "green",
        xmin = wt_temp_min,
        xmax = wt_temp_max,
        ymin = wt_deriv_min,
        ymax = wt_deriv_max
      ) +
      theme(text = element_text(size = 18))

    # Return the plot object
    ggp

  })

  # Calling genotypes
  output$genotypeplot <- renderPlot({

    # Access reactive data
    data <- data()

    # access reactive params
    mut_temp_min <- input$mut_temp_min
    mut_temp_max <- input$mut_temp_max

    wt_temp_max <- input$wt_temp_max
    wt_temp_min <- input$wt_temp_min

    wt_deriv_max <- input$wt_deriv_max
    wt_deriv_min <- input$wt_deriv_min

    mut_deriv_max <- input$mut_deriv_max
    mut_deriv_min <- input$mut_deriv_min


    ## define which samples has peaks in these areas
    peak_mut_positions <-
      data %>%
      group_by(position) %>%
      dplyr::filter(deriv %>% between(mut_deriv_min, mut_deriv_max)) %>%
      dplyr::filter(temp %>% between(mut_temp_min, mut_temp_max)) %>%
      dplyr::filter(lag(deriv, 1) < deriv & lead(deriv, 1) < deriv) %>%
      select(temp, deriv) %>%
      .$position %>% unique

    peak_wt_positions <-
      data %>%
      group_by(position) %>%
      dplyr::filter(deriv %>% between(wt_deriv_min, wt_deriv_max)) %>%
      dplyr::filter(temp %>% between(wt_temp_min , wt_temp_max)) %>%
      dplyr::filter(lag(deriv, 1) < deriv & lead(deriv, 1) < deriv) %>%
      select(temp, deriv) %>%
      .$position %>% unique

    # make a genotype data frame
    genotypes <- tibble(
      position = data$position %>% unique,
      wt_allele = case_when(
        position %in% peak_wt_positions ~ TRUE,
        TRUE ~ FALSE),
      mut_allele = case_when(
        position %in% peak_mut_positions ~ TRUE,
        TRUE ~ FALSE),
      genotype = case_when(
        wt_allele == FALSE & mut_allele == FALSE ~ "NA",
        wt_allele == TRUE & mut_allele == FALSE ~ "wt",
        wt_allele == TRUE & mut_allele == TRUE ~ "het",
        wt_allele == FALSE & mut_allele == TRUE ~ "hom",
      ) %>%
        factor(levels = c("wt", "het", "hom", "NA"))
    )

    # plot
    data %>%
      left_join(genotypes) %>%
      ggplot(
        aes(x = temp, y = deriv, colour = genotype)
      ) +
      geom_line(
        aes(group = position)
      ) +
      scale_y_continuous(labels = comma) +
      scale_x_continuous( # zoom into region of interest.
        limits = c(78, 87),
        breaks = 1:90
      ) +
      scale_color_manual(
        values = c("darkgreen", "red", "darkred","grey50")
      ) +
      facet_wrap(~genotype, nrow = 1) +
      theme(text = element_text(size = 18))

    })

  output$downloadGeno <- downloadHandler(



    filename = "genotypes.csv",
    content = function(file) {

      # Access reactive data
      data <- data()

      # access reactive params
      mut_temp_min <- input$mut_temp_min
      mut_temp_max <- input$mut_temp_max

      wt_temp_max <- input$wt_temp_max
      wt_temp_min <- input$wt_temp_min

      wt_deriv_max <- input$wt_deriv_max
      wt_deriv_min <- input$wt_deriv_min

      mut_deriv_max <- input$mut_deriv_max
      mut_deriv_min <- input$mut_deriv_min


      ## define which samples has peaks in these areas
      peak_mut_positions <-
        data %>%
        group_by(position) %>%
        dplyr::filter(deriv %>% between(mut_deriv_min, mut_deriv_max)) %>%
        dplyr::filter(temp %>% between(mut_temp_min, mut_temp_max)) %>%
        dplyr::filter(lag(deriv, 1) < deriv & lead(deriv, 1) < deriv) %>%
        select(temp, deriv) %>%
        .$position %>% unique

      peak_wt_positions <-
        data %>%
        group_by(position) %>%
        dplyr::filter(deriv %>% between(wt_deriv_min, wt_deriv_max)) %>%
        dplyr::filter(temp %>% between(wt_temp_min , wt_temp_max)) %>%
        dplyr::filter(lag(deriv, 1) < deriv & lead(deriv, 1) < deriv) %>%
        select(temp, deriv) %>%
        .$position %>% unique

      # make a genotype data frame
      genotypes <- tibble(
        position = data$position %>% unique,
        wt_allele = case_when(
          position %in% peak_wt_positions ~ TRUE,
          TRUE ~ FALSE),
        mut_allele = case_when(
          position %in% peak_mut_positions ~ TRUE,
          TRUE ~ FALSE),
        genotype = case_when(
          wt_allele == FALSE & mut_allele == FALSE ~ "NA",
          wt_allele == TRUE & mut_allele == FALSE ~ "wt",
          wt_allele == TRUE & mut_allele == TRUE ~ "het",
          wt_allele == FALSE & mut_allele == TRUE ~ "hom",
        ) %>%
          factor(levels = c("wt", "het", "hom", "NA"))
      )

      write.csv(genotypes, file)
    }
  )



}

# Run the application
shinyApp(ui = ui, server = server)
