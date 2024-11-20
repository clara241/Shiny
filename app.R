#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(bslib)
library(tidyverse)

counts_liver <- read.table("./data/normalizedcounts.tsv", header = T, row.names = 1)
traits <- read.table("./data/traits.txt", header = TRUE, row.names = 1)


pca_data<-prcomp(t(counts_liver), scale=T)
summary(pca_data)
pca_sig.var<-pca_data$sdev^2
pca_sig.var.per<-round(pca_sig.var/sum(pca_sig.var)*100, 1)

pca_sig.data<-data.frame(Sample=rownames(pca_data$x), PC1=pca_data$x[,1], PC2=pca_data$x[,2], PC3=pca_data$x[,3], PC4=pca_data$x[,4], PC5=pca_data$x[,5])
pca_sig.data<-pca_sig.data[-1]




data <- read.table("./data/diffex_out.csv", header = TRUE, sep = ",")
# functions
pvalue_candidate_f <- function(x) {
  if (class(data[[x]]) == "numeric") {
    if (max(data[[x]], na.rm = TRUE) <= 1) {
      if (min(data[[x]], na.rm = TRUE) >= 0) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

gene_candidate_f <- function(x) {
  if (class(data[[x]]) == "character") {
    return(TRUE)}
  return(FALSE)
}

# check data columns to ID most likely candidates for pval, logfc, and gene ID
pval_cols <- names(data)[sapply(names(data), pvalue_candidate_f)]
logfc_cols <- names(data[,2:5])
gene_cols <- names(data)[sapply(names(data), gene_candidate_f)]



# Define UI ----
ui <- page_sidebar(
  title = "OBECAN",
  sidebar = sidebar(
    helpText(
      "Display OBECAN results from bulk-RNA analysis of liver and adipose tissue."
    ),
    
    # Select comparison
    selectInput(inputId = "comparison",
      label = "Choose a comparison to display",
      choices = colnames(traits),
      selected = "Fibrosis"
    ),
    
    # Select number of samples to display in PCA
    sliderInput(inputId = "n_samp", label = "Number of samples to display",
                min = 5, max = 91, value = 10)
    ),

  accordion(  
    accordion_panel( 
      title = "PCA plot", 
      icon = bsicons::bs_icon("menu-app"),
      sidebarPanel(
        selectInput('xcol', 'X Variable', names(pca_sig.data[,1:5])),
        selectInput('ycol', 'Y Variable', names(pca_sig.data[,1:5]))
      ),
      plotOutput(outputId = "plot1", height = "400px", width = "600px")
    ),  
    accordion_panel(
      title = "Volcano plot",
      icon = bsicons::bs_icon("sliders"),
      sidebarPanel(
        # VOLCANO PLOT PANEL -----
        tabPanel("Volcano Plot",
                 h2("Interactive Volcano Plot"),
                 sidebarLayout(
                   
                   # VOLCANO PLOT SIDE PANEL ------
                   sidebarPanel(width = 10,
                                
                                # SELECT AXES LABELS -----
                                h4("Select volcano plot axes:"),
                                
                                # select column for pval
                                selectInput("pvalue_col",
                                            "Input column for significance (y axis)",
                                            pval_cols,
                                            multiple = FALSE),
                                
                                # select column for fold change
                                selectInput("logfc_col",
                                            "Input column for effect size (x axis)",
                                            logfc_cols,
                                            multiple = FALSE),
                                # SET PVAL AND LOGFC THRESHOLDS ----- 
                                h4("Set significance and effect size thresholds:"),
                                
                                # set pvalue threshold 
                                sliderInput("pvalue_threshold",
                                            "Set significance threshold",
                                            min = 0,
                                            max = 1,
                                            value = .05),
                                
                                # set logfc threshold
                                uiOutput("logfc_slider"),
                                
                                
                                # CUSTOMIZE PLOT -----
                                h4("Customize plot:"),
                                
                                # show/hide logfc and pval line
                                checkboxInput("show_pvalue_threshold",
                                              "Show significance threshold line",
                                              value = TRUE),
                                
                                # show/hide logfc lines
                                checkboxInput("show_logfc_threshold",
                                              "Show effect size threshold line",
                                              value = TRUE),
                                
                                # color differentially expressed genes
                                checkboxInput("color_by_de",
                                              "Color significantly different features",
                                              TRUE),
                                
                                # output ui for axis label inputs
                                uiOutput("y_axis_labeler"),
                                uiOutput("x_axis_labeler"),
                                
                                # label legend
                                textInput("legend_title",
                                          "Specify legend title",
                                          value = "Differentially Expressed")),
                   
                   # VOLCANO PLOT MAIN PANEL -----
                   mainPanel(
                     # output info from click
                     p(strong("Plot interactivity:")),
                     p("- View a point's feature label, effect size, and significance by hovering over a point."),
                     p("- Add labels to a feature of interest by using the gene selection dropdown in the sidebar or clicking on the point on the plot."),
                     p("- Remove a label by deleting the selection from the gene selection dropdown or clicking the point on plot a second time."),
                     p("- To zoom click and drag over the plot to select the area you wish to zoom in on. Then, double click to zoom into the selected area. Double click again to zoom out."),
                     verbatimTextOutput("click_info",
                                        placeholder = TRUE),
                     
                     # output ggplot volcano
                     plotOutput("volcano_plot",
                                width = "100%",
                                height = "600px",
                                hover = "volcano_hover",
                                click = "volcano_click",
                                dblclick = "volcano_dbl_click",
                                brush = brushOpts(
                                  id = "volcano_brush",
                                  resetOnNew = TRUE)),
                     
                     # Download button for plot
                     downloadButton('download_volcano', 'Download volcano plot as PDF'),
                     
                     br(),
                     br(),
                     
                     # HIGHLIGHTED GENES TABLE -----
                     dataTableOutput("gene_highlight_tbl"))
                   
                 ) # end sidebarLayout
        ), # end volcano plot tabPanel)
    )),  
    accordion_panel(
      title = "Data",
      icon = bsicons::bs_icon("bar-chart"),
      # DATA PANEL -----
      tabPanel("Data",
               sidebarLayout(
                 
                 # DATA PANEL SIDEBAR
                 sidebarPanel(width = 10,
                              
                              # some text explanation
                              em("Threshold for what is considered differentially expressed is set in Volcano Plot tab by using 
                                 the significance and effect size sliders"),
                              
                              # Show differentiall expressed genes only
                              checkboxInput("show_de",
                                            "Show only significantly different features",
                                            FALSE)),
                 
                 # DATA PANEL MAIN PANEL
                 mainPanel(dataTableOutput("gene_data")))
      ) # end data tab panel
    ),  
    accordion_panel(
      title = "Section D",
      icon = bsicons::bs_icon("calendar-date"), 
      "Section D content" 
    ),  
    id = "acc",  
    open = NULL  
  )
  

)



# Define server logic ----
server<- function(input, output, session) {
  # Reactive expression to update 'category' based on selected comparison
  selectedData <- reactive({
    # Update the 'category' column based on the selected comparison
    pca_sig.data2 <- pca_sig.data
    pca_sig.data2$category <- traits[[input$comparison]]
    pca_sig.data2[c(input$xcol, input$ycol, "category")]
  })
  
  # Create new df that is n_samp obs from selected type movies
  data_filtered <- reactive({
    req(input$n_samp)
    sample_n(selectedData(), input$n_samp)
  })
  
  output$plot1 <- renderPlot({
    
    ggplot(data=data_filtered(), aes(x=.data[[input$xcol]], y = .data[[input$ycol]], colour=category))+
      geom_point(size=2, stroke=1, alpha=0.8, aes(color=category))+
      xlab(input$xcol) +
      ylab(input$ycol) +
      theme_bw()+
      ggtitle("My PCA Graph")
  })
  
  
 
  
  
  
  
}



# Run the app ----
shinyApp(ui = ui, server = server)


