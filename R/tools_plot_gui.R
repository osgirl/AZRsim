###############################################################################
###############################################################################
# This file contains at the moment a simple shiny app to display the contents
# of a dataframe that somes back from the simulation of an AZRmodel.
###############################################################################
###############################################################################

###############################################################################
# AZRplot: Use shiny to graphically represent the contents of a dataframe
###############################################################################
#' Use shiny to graphically represent the contents of a dataframe
#'
#' As input data frames are used and the columns plotted.
#' By default it looks for a "time" to be used as X-axis.
#' In the future it might handle "ID" columns in some way and others,
#' such as "TAD" and "AMT", ... Lets start simple for now.
#' Only numeric columns in the dataframe are considered.
#'
#' @param data A dataframe
#' @export

shiny_plot <- function(data) {

  # Check dataframe
  if (!is.data.frame(data))
    stop("AZRplot: 'data' is not a data frame")

  # Keep "numeric" columns only
  dataNum <- data[,unname(which(sapply(data, class)=="numeric"))]

  # Check if a TIME column is present if yes then move it to first column
  # Use case insensitive matching
  indexTIME <- unname(which(toupper(names(dataNum))=="TIME"))

  # Reorder columns in dataframe if "TIME" was found
  if (length(indexTIME)>0)
    dataNum <- dataNum[,c(indexTIME,setdiff(1:ncol(dataNum),indexTIME))]

  #################################################################
  # Create the shiny app
  # dataNum: has to be a data frame with numeric columns only.
  #          First column will be used as default x-axis
  #################################################################

  shiny::shinyApp(

    #########################################
    # shiny ui function
    #########################################

    ui = shiny::fluidPage(
      shiny::titlePanel("AZRplot (simple exploration)"),

      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("Xname","X-Axis",choices = names(dataNum),selected=names(dataNum)[1]),
          shiny::checkboxGroupInput("Yname","Y-Axis",choices = names(dataNum),selected=names(dataNum)[2])
        ),
        shiny::mainPanel(shiny::plotOutput("mainPlot"))
      )
    ),

    #########################################
    # shiny server function
    #########################################

    server = function(input, output) {
      output$mainPlot <- shiny::renderPlot({
        # Keep only selected Ynames
        dat <- dataNum[,unique(c(input$Xname,input$Yname))]
        if (length(names(dat)>1)) {
          dat <- tidyr::gather(dat,"variable","value",-eval(parse(text=input$Xname)))
          p <- ggplot2::ggplot(data=dat,ggplot2::aes(x=dat[,input$Xname],y=dat[,"value"],group=dat[,"variable"],color=dat[,"variable"])) +
            ggplot2::geom_line(size=1) +
            ggplot2::labs(x = input$Xname) +
            ggplot2::labs(y = NULL) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text=ggplot2::element_text(size=12),axis.title=ggplot2::element_text(size=14,face="bold"))
          p
          return(p)
        } else {
          return(NULL)
        }
      })
    }
  )
}
