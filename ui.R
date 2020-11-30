library(shiny)
library(Seurat)
library(shinythemes)
library(shinycustomloader)
shinyUI(fluidPage(
  theme=shinytheme("cerulean"),
  #themeSelector(),

  
  navbarPage(title="Cancer Smell",
                   tabPanel("ABOUT",
                          h4("A computational workflow that systematically estimates the activation status of the chemosensory olfactory
and taste receptors at the single-cell resolution."),
                            HTML('<center><img src="RStudio-Ball_2.png", heigth=400, width=400></center>')
                            #tags$img(src='RStudio-Ball_2.png', heigth=400, width=400,)
                            
                            ),
                   tabPanel("DATA",
                            sidebarLayout(
                              sidebarPanel(
                                fileInput("file","Upload the file"),
                                                        h4("Select Parameters"),
                                                          checkboxInput(inputId = 'header', label = 'Header', value = FALSE),
                                                          checkboxInput(inputId = "stringAsFactors", "stringAsFactors", FALSE),
                                                          radioButtons(inputId = "rowname", label = "Row Names", choices=c(Yes=1, No=NULL), 1),
                                                          radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ',')
                                                          #radioButtons("choice","Choose an option", choices=c("Dataset" = 1, "Structure" = 2,"Summary" = 3 ))
                                       
                                
                              ),
                              mainPanel(tableOutput("dat"),
                                        verbatimTextOutput("dimen")
                            )
                            
                            )),
             tabPanel("OBJECT",
                      actionButton("action0","Create Object"),
                      withLoader(verbatimTextOutput("seu_object"))
                      
                      ),
             tabPanel("DOWNSTREAM ANALYSIS",
                      #h4("Visualize the Feature Plot"),
                      actionButton("action1","Visualize Variable Features"),
                      withLoader(plotOutput("p1")),
                      #downloadButton("p1_down",label="Download as PDF"),
                      #submitButton("Next"),
                      #h4("Normalize Data"),
                      actionButton("action2","Normalize Data"),
                      withLoader(verbatimTextOutput("data_norm")),
                      verbatimTextOutput("variable_features"),
                      actionButton("action3","Plot Variable Features"),
                      withLoader(plotOutput("plot_vf")),
                      actionButton("action4","Scale Data"),
                      verbatimTextOutput("scale_data"),
                      actionButton("action5","Run PCA"),
                      withLoader(plotOutput("pca")),
                      actionButton("action6","Plot HeatMap"),
                      withLoader(plotOutput("heatmap")),
                      actionButton("action7","Jack Straw Plot"),
                      withLoader(plotOutput("jackstraw")),
                      actionButton("action8","Elbow Plot"),
                      withLoader(plotOutput("elbow")),
                      actionButton("action9","Find Neighbours"),
                      withLoader(verbatimTextOutput("neigh")),
                      actionButton("action10","Run UMAP"),
                      withLoader(plotOutput("umap")),
                      actionButton("action11","Find Markers"),
                      withLoader(verbatimTextOutput("markers")),
                      actionButton("action12","Display Markers Table"),
                      tableOutput("markers_display")
                      
                      
             ),
             tabPanel("ORs",
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("f_or")
                          
                        ),
                        mainPanel(plotOutput("or_umap"))
                      )
                      #tableOutput("f_or"))
                   ),
             tabPanel("Other Receptors",
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("other_or")
                          
                        ),
                        mainPanel(plotOutput("other_umap"))
                      )
                      #tableOutput("f_or"))
             ),
             
             tabPanel("User Input",
                      sidebarLayout(
                        sidebarPanel(("Enter gene name"),
                          textInput("usergene","","")
                        ),
                        mainPanel(plotOutput("user_umap"))
                      )
                      #tableOutput("f_or"))
             )
)
)
)
