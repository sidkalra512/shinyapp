library(shiny)
library(ggplot2)
library(Seurat)
library(shinythemes)
options(shiny.maxRequestSize = 500*1024^2)

shinyServer(function(input,output){
  
  # This reactive function will take the inputs from UI.R and use them for read.table() to read the data from the file. It returns the dataset in the form of a dataframe.
  # file$datapath -> gives the path of the file
  data <- reactive({
    file1 <- input$file
    if(is.null(file1)){return()} 
    read.table(file=file1$datapath, sep=input$sep, header = input$header, stringsAsFactors = input$stringAsFactors, row.names = as.numeric(input$rowname))
    
  })
  
  output$dat<-renderTable({
    if(is.null(data())){return ()}
    data()[1:10,1:7]
  })
  
  output$dimen <- renderPrint({
    print("Dimensions")
    dim(data())
  })

  
  #pb <- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200)
  
  output$seu_object <- renderPrint({
    if(input$action0==0){
      return()
    }
    pb <<- CreateSeuratObject(counts = data(), min.cells = 3, min.features = 200)
    print(pb)
    print("Your object has been created, start downstream analysis now!!!!")
  })
  
  output$p1<-renderPlot({
    if(input$action1==0){
      return()
    }
      
    
    VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    
  })
  
  output$p1_down<-downloadHandler(
    filename<<-"Featureplot.pdf",
    content=function(filename){
      pdf(file=filename)
      VlnPlot(pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
      dev.off()
    }
  )
  
  output$data_norm<-renderPrint({
    if(input$action2==0){
      return()
    }
    pb <<- NormalizeData(pb)
    print("Normalization Done")
  }
  )
  
  
  output$variable_features<-renderPrint({
      pb <<- FindVariableFeatures(pb, selection.method = "vst", nfeatures = 2000)
      #list_of_variable_features<<-VariableFeatures(pb)
      #list_of_variable_features<<-as.data.frame(list_of_variable_features)
      top10 <<- head(VariableFeatures(pb), 20)
  })
  
  output$plot_vf<-renderPlot({
    if(input$action3==0){
      return()
    }
      plot1 <- VariableFeaturePlot(pb)
      plot2 <- LabelPoints(plot = plot1, points = top10,repel=TRUE,xnudge = 0,ynudge = 0)
      plot(plot2)
  })

  output$scale_data<-renderPrint({
    if(input$action4==0){
      return()
    }
    all.genes <- rownames(pb)
    pb <<- ScaleData(pb, features = all.genes)
    print("Centering and scaling data matrix")
    print("100% Done")
  })
    
    output$pca<-renderPlot({
      if(input$action5==0){
        return()
      }
      pb <<- RunPCA(pb, features = VariableFeatures(object = pb))
      DimPlot(pb, reduction = "pca")
    })
    
    output$heatmap<-renderPlot({
      if(input$action6==0){
        return()
      }
      DimHeatmap(pb, dims = 1:6, cells = 500, balanced = TRUE)
      
    })
    
    output$jackstraw<-renderPlot({
      if(input$action7==0){
        return()
      }
      pb <<- JackStraw(pb, num.replicate = 100)
      pb <<- ScoreJackStraw(pb, dims = 1:20)
      #pdf(file="Main_pipeline/GSE756881/Jackstraw_plot.pdf")
      JackStrawPlot(pb, dims = 1:15)
    })
    
    output$elbow<-renderPlot({
      if(input$action8==0){
        return()
      }
      ElbowPlot(pb)
      
    })
    
    output$neigh<-renderPrint({
      if(input$action9==0){
        return()
      }
      pb <<- FindNeighbors(pb, dims = 1:10)
      pb <<- FindClusters(pb, resolution = 0.5)
      head(Idents(pb), 5)
    })
    
    output$umap<-renderPlot({
      if(input$action10==0){
        return()
      }
      pb <<- RunUMAP(pb, dims = 1:10)
      DimPlot(pb, reduction = "umap")
      
    })
    
    output$markers<-renderPrint({
      if(input$action11==0){
        return()
      }
        pb.markers <<- FindAllMarkers(pb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      
    })
    
    output$markers_display<-renderTable({
      if(input$action12==0){
        return()
      }
     head(pb.markers)
    })
    
    output$f_or<-renderUI({
      or<-read.csv("functional_or.csv")
      row_names<-row.names(pb)
      row_names<-as.data.frame(row_names)
      colnames(row_names)<-"Receptor_names"
      factor0<<-merge(row_names,or,by.x="Receptor_names",by.y="Symbol")
      factor0
      selectInput("detected_or","Select Receptor",choices = factor0)
      #paste("Number of ORs detected",nrow(factor0))
    })
    
    output$or_umap<-renderPlot({
      FeaturePlot(pb, features = c(input$detected_or))
    })
    
    output$other_or<-renderUI({
      other_or<-read.csv("other_receptors.csv")
      other_receptors<-row.names(pb)
      other_receptors<-as.data.frame(other_receptors)
      colnames(other_receptors)<-"Receptor_names"
      factor_other<<-merge(other_receptors,other_or,by.x="Receptor_names",by.y="Symbol")
      selectInput("detected_other_or","Select Receptor",choices = factor_other)
      #paste("Number of ORs detected",nrow(factor0))
    })
    
    output$other_umap<-renderPlot({
      FeaturePlot(pb, features = c(input$detected_other_or))
    })
    
    output$user_umap<-renderPlot({
      FeaturePlot(pb, features = c(input$usergene))
    })

})