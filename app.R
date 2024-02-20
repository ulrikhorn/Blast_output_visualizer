
library(shiny)
library(tidyverse)
library(magrittr)
library(rentrez)
library(gtools)
library(plotly)
library(htmlwidgets)
#library(sys)
#Sys.setlocale("LC_ALL", "en_US.UTF-8")

ui <- fluidPage(

    # Application title
    titlePanel("Visualization of Blast output"),
    

    # Sidebar with inputs
    sidebarLayout(
        sidebarPanel(
          fileInput("file", "Select the blastn output file"),
          fileInput("meta", "Make a file of the blastn --outputfmt argument used, in a text file and upload here"),
          uiOutput("contig_select"),
          numericInput("select_bitscore", "Select Bitscore Threshold", min = 0, max = 50000, value = 10000),
          downloadButton("report", "Generate report")),

        # Show a plot of the generated distribution
        mainPanel(
          plotlyOutput("contig_plot"),
          plotlyOutput("pie_shart"),
          plotlyOutput("length_plot"),
          plotOutput("length_hist"),
          plotlyOutput("bitscore_plot")
        )
    )
)

# Define server logic required 
server <- function(input, output) {

  ## Data wrangeling for contig alignment visualization
  data <- reactive({
    req(input$file)
    req(input$meta)
    raw_df <- read_tsv(input$file$datapath, col_names = str_split(names(read_tsv(input$meta$datapath)), pattern = " ")[[1]][-1])
    raw_df$acc_nr <- raw_df$sseqid %>% str_replace("\\.1_.*", "")
    
    org_list <- vector()
    for (i in unique(raw_df$acc_nr)) {
      org_list <-
        append(org_list, entrez_summary(db = "nucleotide", id = i)$organism)
    }
    ncbi_df <-
      tibble(acc_nr = unique(raw_df$acc_nr), organism = org_list)
    raw_df <- left_join(raw_df, ncbi_df)
    
    
    raw_df <- distinct(raw_df)
    
    raw_df$nr <- seq(from = 1, to = nrow(raw_df))
    raw_df <- raw_df %>% group_by(organism) %>% mutate(no_matches = n())
    raw_df
    })
  
  # Creating a list of contig options from the dataset
  output$contig_select <- renderUI({
    contig_data <- data()
    contig_list <- mixedorder(unique(contig_data$qseqid))
    selectInput("contig_list", "Select contig", contig_list)
  })
  
  # Creating an object that gets the selected contig
  contig_selection <- reactive({
    req("contig_list")
    selected_contig <- input$contig_list
    selected_contig
  })
  
  #Creating an object to retreive the selected bitscore threshold
  bitscore_selection <- reactive({
    req(input$select_bitscore)
    selected_bitscore <- as.numeric(input$select_bitscore)
    selected_bitscore
  })
  
  #Ploting the contig alignment
    output$contig_plot <- renderPlotly({
      ggplotly(ggplot(data() %>% filter(qseqid == paste("contig_", contig_selection(), sep ="") & bitscore > bitscore_selection()), aes(x = qstart, xend= qend, y = as.character(nr), yend= as.character(nr), text = organism, color = bitscore))+
        geom_segment(position = "stack", size = 3)+
        geom_segment(mapping = aes(x = 0, xend = qlen, y = qseqid, yend = qseqid), size = 3, color = "red")+
          theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank())+
        guides(fill = FALSE), tooltip = c("qseqid", "organism", "qstart", "qend", "bitscore", "nr")) %>% layout(hovermode = "x")
        
    })
    
    ## Pie chart visualization
    pie_data <- reactive({
      req(input$file)
      req(input$meta)
      pie_df <- data() %>% mutate(percent_matches = no_matches/nrow(data()))
      pie_df
    })
    
    output$pie_shart <- renderPlotly({
      plot_ly(data = pie_data(), values = ~no_matches, labels=~organism, type = "pie")
    })
    
    
    length_data <- reactive({
      req(input$file)
      req(input$meta)
      length_df <- data() %>% group_by(qseqid, qlen) %>% summarise()
      length_df <- length_df[order(length_df$qlen),]
      level = rev(length_df$qseqid)
      length_df
    })
    
    level_data <- reactive({
      req(input$file)
      req(input$meta)
      level_df <- length_data()
      level = rev(level_df$qseqid)
      level
    })
    
    output$length_plot <- renderPlotly({
      ggplotly(ggplot(length_data())+
        geom_col(mapping = aes(factor(qseqid, levels = level_data()), qlen, fill = qseqid))+
        theme(axis.text.x = element_text(angle = 90))+
        guides(fill = FALSE)+
        labs(x = "Contig", y = "Contig Length (bp)"), tooltip = c("qseqid", "qlen"))
    })
    
    output$bitscore_plot <- renderPlotly({
      ggplotly(ggplot(data())+
                 geom_col(mapping = aes(qseqid, median(bitscore)))+
                 theme(axis.text.x = element_text(angle = 90)), tooltip = "qseqid") %>% layout(hovermode = "x") 
    })
    
    output$length_hist <- renderPlot({
      ggplot(data() %>% filter(qseqid == "contig_3"), aes(length, fill = nr))+
        geom_histogram()
    })
    
    output$report <- downloadHandler(
      filename = "blast_report.html",
      content = function(file){
        params <- list(blast_output = input$file$datapath, blast_outputformat = input$meta$datapath, shiny_data = data(), contig = contig_selection(), bitscore = bitscore_selection())
        rmarkdown::render("blast_report.Rmd", 
                          output_file = file,
                          params = params
                          #envir = new.env(parent = globalenv())
                          )
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
