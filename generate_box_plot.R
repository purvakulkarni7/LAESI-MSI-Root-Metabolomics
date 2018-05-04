#' generate_box_plot.R
#' This function takes a csv file as input that contains normalized values of the sample concentrations for the replicates of a pair of plant species  (ex. CS and CJ)
#'                                                      
#' @param filePath path to the *.csv file that contains the normalized values of the sample concentrations. This file can be exported from metaboAnalyst (http://metaboanalyst.ca) ver 4.0 after performing normalization. 
#' 
#' @param mass the mass value for which the box plot has to be generated. 
#'
#' @export Exports *.png file with a box plot for that specific mass. 
#' 
#' @note The files and masses used as an input for this code are:
#' metaboanalystdata_normalizedGMGP.csv
#' "158.2646933"
#' "172.3828554"
#' "196.5855371"
#' "250.8271064"
#' 
#' @author Purva Kulkarni
#' 

generate_box_plot <- function(filePath, mass)
{
  #Read file data
  boxplotData <- read.table(filePath, sep = ",")
  
  # Remove first column from the data frame
  boxplotData[,1] <- NULL
  
  # Assign first row as column names and ttranspose the data.frame using the given steps
  colnames(boxplotData) <- boxplotData[1,]
  boxplotData = boxplotData[-1, ]
  n <- boxplotData$`3`
  n <- as.character(n)
  boxplotData <- as.data.frame(t(boxplotData[,-1]))
  colnames(boxplotData) <- n
  
  # Subset this data.frame to extract the row with the input mass value
  boxplotDataSubset <- subset(boxplotData,  row.names.data.frame(boxplotData) == mass)
  
  # Arrange values for two species in individual columns
  col1 <- rbind(boxplotDataSubset[,1], boxplotDataSubset[,2], boxplotDataSubset[,3])
  col2 <- rbind(boxplotDataSubset[,4], boxplotDataSubset[,5], boxplotDataSubset[,6])
  plotData <- cbind(col1,col2)
  colnames(plotData) <- c(colnames(boxplotData[1]),colnames(boxplotData[4]))
  
  # Generate the boxplot
  par(bg=NA)
  boxplot(plotData)
  dev.copy(png,'boxplot_plot.png')
  dev.off()
}