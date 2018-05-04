#' exportMassSpectrumObjects.R
#' This function reads all the mass spectrum text files present in a folder and generates a list of mass spectum objects from it. The length of this list is equal to the number of mass spectrum files present in the folder.  
#'                                                      
#' @param inputPath input folder path containing text files. Each text file contain only two column with mass and intensity values
#' 
#' @return An object of class list. This object contains a list of mass spectrum object (An object type from the package maldiquant)
#' 
#' @note This script requires two packages MALDIquant, MALDIquantForeign
#'
#' @author Purva Kulkarni
#'

library(MALDIquant)
library(MALDIquantForeign)

massFeatures <- function(inputPath)
{
  fileList = list.files(path = inputPath, pattern = "*.txt", full.names = FALSE)
  setwd(inputPath)
  
  massSpectrumObjects <- vector(mode="numeric", length=length(fileList))
  
  # Read all the files one by one
  for (i in 1:length(fileList))
  {
    massSpectrumObjects[i] <- import(fileList[i])
  }
  
  return(massSpectrumObjects)
}