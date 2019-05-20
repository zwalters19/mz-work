#Zach Walters
#App for visualization of MS-FInder output

library(plotly)
library(shiny)
library(clipr)
library(stringr)
library(matrixTests)
library(reticulate)


server = function(input, output, session){
  options(shiny.maxRequestSize=6*1024^2) 
  options(rsconnect.max.bundle.size= 4*1024^2)
  options(shiny.reactlog=TRUE)
  options(shiny.loadFactor = .9)
  
  #a trio of functions that perform more cpu-intensive calculations. Invoked by data_intermediate, should 
  #only be called when input files change
  data_parsing <- function(files){
    count = 0
    for(i in (files)){
      if(count == 0){
        data_table <- read.csv(i, sep = "\t", header = T,stringsAsFactors = F)
        intensity <- as.numeric(data_table[,16])
        count = count + 1
      }
      else{
        data_table <- read.csv(i, sep = "\t", header = T,stringsAsFactors = F)
        intensity = cbind(intensity, as.numeric(data_table[,16]))
        
      }
    }
    return(intensity)
  }
  t_test_fun <- function(int1, int2){
    if(is.matrix(int1) & is.matrix(int2)){
      results <- row_t_welch(int1,int2)
    }
    return(results)
  }
  dysregulated <- function(int1, int2){
    max1 <- apply(int1,1,mean)
    max2 <- apply(int2,1,mean)
    updown <- max1/max2
    return(updown)
  }
  
  #first-pass data processing. Since this is the most computationally-expensive portion, 
  #should only be called upon file input
  data_intermediate <- reactive({
    out_temp <-  read.csv(input$file1$datapath[[1]], sep = "\t", header = T,stringsAsFactors = F)
    out <-  out_temp[,-ncol(out_temp)]
    int1 <- data_parsing(input$file1$datapath)
    out$Intensity <- apply(int1,1,max)
    if(!is.null(input$file2)){
      int2 <- data_parsing(input$file2$datapath)
      int_comb <- cbind(int1, int2)
      out$Intensity <- apply(int_comb, 1, max)
      out$IntMean <- apply(int_comb, 1, mean)
      out$fold <- dysregulated(int1, int2)
      out$updown <- out$fold > 1
      statistics_res <- t_test_fun(int1,int2)
      out$p_value <- statistics_res$pvalue
      out$t_value <- statistics_res$statistic
    }
    else{
      out$fold <- rep(0,nrow(out))
      out$p_value <- rep(0,nrow(out))
      out$updown <- rep(F,nrow(out))
    }
    
    out$mscategory  <-  out$ms2_match != ""
    out$ms2id[out$MS2Match != ""] <- out$MS2Name[out$MS2Match != ""]
    out$ms2id[out$MS2Match == ""] <- out$Formula[out$MS2Match == ""]
    if(input$box1 == T){
      formulas = unique(out$Formula)
      for(f in formulas){
        matching = subset(out, out$Formula) == f
        if(nrow(matching)>1){
          matching <- matching[order(matching$Intensity),]
          out <- out[!matching[-1,]]
        }
        
      }
    }
    out$compound_score[is.na(out$compound_score)] <- 0
    category1 <- out$formula == ""
    category2 <- out$inchi == ""
    out$category <- category1 + category2
    return(out)
  })
  
  #Subsetting based on parameters.Updates every time parameters change, but fast
  data <- reactive({
    out <- data_intermediate()
    if(!is.null(input$file2)){
      out <- subset(out,out$fold > input$text2 | out$fold < 1/input$text2)
      out <- subset(out, out$p_value < input$text3)
      
    }
    out <- subset(out, out$Intensity > input$text1)
    out <- subset(out, out$ms1_score > input$text4)
    
    if(input$formula_name != ""){
      out = subset(out, out$formula == input$formula_name)
    }
    return(out)
  })
  ms2_data <- reactive({
    out <- data()
    out <- subset(out, out$ms2_source != "")
    out <- subset(out, out$compound_score > input$text5)
    
    return(out)
  })
  
  
  #plot ms1 data 
  output$plot <- renderPlotly({
    req(input$file1)
    out <- data()
    y = list(
      title = "m/z"
    )
    x = list(
      title = "RT"
    )
    
    plot_ly(out, type = "scatter",x = out$rt, y = out$mz, 
            text =  out$formula,mode = "markers", marker = list(size = 3, width = 3),
            color = out$updown, hoverinfo = "text", source = "plot1",
            colors=c('red','black','darkgreen'))%>%
      layout(yaxis = y, xaxis = x, showlegend = F, title = "MS1")%>%
      hide_colorbar()
    
  })
  output$plot2 <- renderPlotly({
    req(input$file1)
    out <- ms2_data()
    y = list(
      title = "m/z"
    )
    x = list(
      title = "RT"
    )
    plot_ly(out, type = "scatter",x = out$rt, y = out$mz,  
            text =  out$compound_name,mode = "markers", color = out$mscategory,source = "plot2",split = out$mscategory,
            hoverinfo = 'text',colors = colorRamp(c("gold", "blue"))(0:2/2))%>%
      layout(yaxis = y, xaxis = x, showlegend = F, title = "MS2")%>%
      hide_colorbar()
  })
  output$click1 <- renderPrint({
    req(input$file1)
    c <- event_data(event = "plotly_click", source = "plot1")
    if(!is.null(c)){
      out = data()
      up <- subset(out, out$updown == T)
      down <- subset(out, out$updown == F)
      
      if((nrow(down)!=0 & nrow(up) != 0)){
        if(c[[1]] == 0){
          out = down
        }
        else{
          out = up
        }
      }
      output$mz <- renderText({
        text = paste("mz: ",out$ms1_score[c[[2]]+1])
        return(text)
      })
      output$rt <- renderText({
        text = paste("rt: ",out$compound_score[c[[2]]+1])
        return(text)
      })
      output$fold <- renderText({
        text = paste("fold: ",out$fold[c[[2]]+1])
        return(text)
      })
      output$pv <- renderText({
        text = paste("p: ",out$p_value[c[[2]]+1])
        return(text)
      })
      output$maxint <- renderText({
        text = paste("max int: ",out$Intensity[c[[2]]+1])
        return(text)
      })
      output$meanint <- renderText({
        text = paste("mean int: ",out$IntMean[c[[2]]+1])
        return(text)
      })
      y = list(
        title = "intensity"
      )
      x = list(
        title = "isotope"
      )
      write_clip(out$formula[c[[2]]+1])
      iso_source = out$isotope_source[c[[2]]+1]
      iso_match = out$isotope_match[c[[2]]+1]
      text1 = str_split(iso_source,";",simplify = T)
      text2 = str_split(iso_match,";",simplify = T)
      intensity1 = as.numeric(text1)
      intensity2 = as.numeric(text2)
      p3 <- plot_ly(out, type = "scatter", x = c(), y = c(),mode = "markers",colors = c(1))%>%
        layout(yaxis = y, xaxis = x, showlegend = F, title = paste(c(out$formula[c[[2]]+1],
                                                                     paste(c("parent Peak:", out$mz[c[[2]]+1], "m/z"),collapse = " "),
                                                                     paste('Adduct:',out$adduct[c[[2]]+1], sep = " ")), collapse = "\n"))%>%
        hide_colorbar()
      p3 <- add_segments(p3,x= 3,xend = 3, y = 0, yend = 0)
      p3 <- add_segments(p3,x= -1,xend = -1, y = 0, yend = 0)
      output$plot3 <- renderPlotly(p3)
      for( i in 1:3){
        p3 <- add_segments(p3, x = i-1 ,xend = i-1, y = 0, yend = intensity1[i],
                           line = list(color = 'black', width = 2))
        p3 <- add_segments(p3, x = i-1 ,xend = i-1, y = 0, yend = -1*intensity2[i],
                           line = list(color = 'red', width = 1))
      }
      #p = out$Formula[c[[2]]+1]
      #return(p)
      output$table1 <- renderTable({
        cpds = str_replace(out$possible_compounds[c[[2]]+1], coll("["),"")
        cpds = str_replace(cpds, coll("]"),"")
        table = matrix(str_split(cpds,", ",simplify = T),ncol = 1)
        colnames(table) = "Possible Matches"
        return(table)
      })
    }
  })
  
  output$click2 <- renderPrint({
    req(input$file1)
    c <- event_data(event = "plotly_click", source = "plot2")
    if(!is.null(c)){
      out <- ms2_data()
      has_match <- subset(out, out$mscategory == 1)
      no_match <- subset(out, out$mscategory == 0)
     
      if(c[[1]] == 0){
        if(nrow(no_match) > 0){
          z = 0
        }
        else{
          z = 1
        }
      }
      if(c[[1]] == 1){
        z = 1
      }
      y = list(
        title = "intensity"
      )
      x = list(
        title = "m/z",
        range = c(1,1000)
      )
    
      if(z == 0){
        out <- no_match
        output$mz <- renderText({
          text = paste("mz: ",out$ms1_score[c[[2]]+1])
          return(text)
        })
        output$rt <- renderText({
          text = paste("rt: ",out$compound_score[c[[2]]+1])
          return(text)
        })
        output$fold <- renderText({
          text = paste("fold: ",out$fold[c[[2]]+1])
          return(text)
        })
        output$pv <- renderText({
          text = paste("p: ",out$p_value[c[[2]]+1])
          return(text)
        })
        output$maxint <- renderText({
          text = paste("max int: ",out$Intensity[c[[2]]+1])
          return(text)
        })
        output$meanint <- renderText({
          text = paste("mean int: ",out$IntMean[c[[2]]+1])
          return(text)
        })
        text <- str_split(str_split(text,";",simplify = T),",",simplify=T)
        p3 <- plot_ly(out, type = "scatter", x = c(), y = c(),mode = "markers",colors = c(1))%>%
          layout(yaxis = y, xaxis = x, showlegend = F, title = paste(no_match$formula[c[[2]]+1],paste(c("parent Peak:", no_match$mz[c[[2]]+1], "m/z"),collapse = " "), sep = "\n"))%>%
          hide_colorbar()
          for(i in(1:length(text[,1]))){
            p3  <- add_segments(p3,x = text[i,1], y = 0, xend = text[i,1],yend = text[i,2],colors = c(1))
          }
        output$plot3 <- renderPlotly(p3)
      }
      if(z == 1){
        out <- has_match
        output$mz <- renderText({
          text = paste("mz: ",out$ms1_score[c[[2]]+1])
          return(text)
        })
        output$rt <- renderText({
          text = paste("rt: ",out$compound_score[c[[2]]+1])
          return(text)
        })
        output$fold <- renderText({
          text = paste("fold: ",out$fold[c[[2]]+1])
          return(text)
        })
        output$pv <- renderText({
          text = paste("p: ",out$p_value[c[[2]]+1])
          return(text)
        })
        output$maxint <- renderText({
          text = paste("max int: ",out$Intensity[c[[2]]+1])
          return(text)
        })
        output$meanint <- renderText({
          text = paste("mean int: ",out$IntMean[c[[2]]+1])
          return(text)
        })
        
        
        text <- has_match$ms2_source[c[[2]]+1]
        text <- str_split(str_split(text,";",simplify = T),",",simplify=T)
        text2 <- has_match$ms2_match[c[[2]]+1]
        text2 <- str_split(str_split(text2,";",simplify = T),",",simplify=T)
        
        intensity = as.numeric(text[,2])
        intensity = intensity / max(intensity)
        intensity = intensity * 100
        
        intensity2 = as.numeric(text2[,2])
        intensity2 = intensity2 / max(intensity2)
        intensity2 = intensity2 *-100
   
        p3 <- plot_ly(out, type = "scatter", x = c(), y = c(),mode = "markers",colors = rep(0,8))%>%
          layout(yaxis = y, xaxis = x, showlegend = F, title = paste(c(has_match$compound_name[c[[2]]+1],paste(c("parent Peak:", has_match$mz[c[[2]]+1], "m/z"),collapse = " "),has_match$formula[c[[2]]+1]), collapse = "\n"))%>%
          hide_colorbar()
        for(i in(1:length(text[,1]))){
          if(intensity[i] > 1)
          p3  <- add_segments(p3,x = text[i,1], y = 0, xend = text[i,1],yend = intensity[i],
                              line = list(color = 'black',width = 2 ))
        }
        if(input$box1 == T){
          for(i in(1:length(text2[,1]))){
            if(intensity2[i] < -1)
              p3  <- add_segments(p3,text = text2[i,3], x = text2[i,1],
                                  y = 0, xend = text2[i,1],yend = intensity2[i],
                                  line = list(color = 'red',width = 1 ))
          }
        }
        else{
          for(i in(1:length(text2[,1]))){
            if(intensity2[i] < -1)
            p3  <- add_segments(p3,x = text2[i,1], y = 0, xend = text2[i,1],yend = intensity2[i],
                                line = list(color = 'red',width = 1 ))
          }
        }
        output$plot3 <- renderPlotly(p3)
      }
    }
    
  })
    
    output$mummichog <- function(){
      req(input$file1)
      req(input$file2)
      use_condaenv('py2')
      out <- data_intermediate()
      filename = "mummichog.txt"
      out_file <- out[c('mz','rt', 'p_value', 't_value')]
      write.table(out_file, filename, col.names = T, row.names = F, sep = "\t", quote = F)
      mummichog_file = paste(getwd(),
                          '\\Dependencies\\Python\\miniconda2\\envs\\py2\\Lib\\site-packages\\mummichog\\main.py',
                          sep = "")
      py_run_file(mummichog_file)
    }
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(width = 3,
      fileInput("file1","Group 1", multiple = T, accept = ".csv"),
      fileInput("file2", "Group 2", multiple = T, accept = ".csv"),
      numericInput("text1", "Intensity Threshold", value = 10000,min = 0, step = 1000),
      numericInput("text2", "Fold Change", value = 1.5, min = 0, step = .5),
      checkboxInput("box1", "Combine Adducts"),
      numericInput("num1", "retention time threshold", value = .5, min = 0, step = .5),
      numericInput("text3", "p value threshold", value = .05,min = 0,max = 1, step = .01),
      numericInput('text4', 'ms1 score threshold', value = 0,min=0,max=5, step = .5),
      numericInput('text5', 'ms2 score threshold',value = 0, min=0,max=10,step=.5),
      actionButton("mummichog", "Run Mummichog"),
      textInput('formula_name',"formula name search", placeholder = "")
    ),
    mainPanel(
      verbatimTextOutput("click1",placeholder = F),
      verbatimTextOutput("click2",placeholder = F),
      verbatimTextOutput('mz',placeholder = T),
      verbatimTextOutput('rt',placeholder = T),
      verbatimTextOutput('fold',placeholder = T),
      verbatimTextOutput('pv',placeholder = T),
      verbatimTextOutput('maxint',placeholder = T),
      verbatimTextOutput('meanint',placeholder = T),
      
      fluidRow(
        tabsetPanel(
          tabPanel(title = 'ms1',plotlyOutput("plot")),
          tabPanel(title = 'ms2',plotlyOutput("plot2"))
        )
      ),
      fluidRow(
        column(12,plotlyOutput("plot3")),
        column(12,tableOutput('table1'))
      )
    )
  )
)


shinyApp(ui, server)

