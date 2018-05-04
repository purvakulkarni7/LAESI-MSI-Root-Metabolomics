#' pairwiseAnalysis-CJCS.R
#' This function contains mass spectrmetry data analysis steps for the pairwise comaprative analysis of plant species CJ and CS. The function parameters used in this function as customized as per the data.
#' 
#' @param CJ_Replicate_1_path folder path containing 50 text files with mass spectra for CJ replicate 1
#' @param CJ_Replicate_2_path folder path containing 50 text files with mass spectra for CJ replicate 2
#' @param CJ_Replicate_3_path folder path containing 50 text files with mass spectra for CJ replicate 3
#' 
#' @param CS_Replicate_1_path folder path containing 50 text files with mass spectra for CS replicate 1
#' @param CS_Replicate_2_path folder path containing 50 text files with mass spectra for CS replicate 2
#' @param CS_Replicate_3_path folder path containing 50 text files with mass spectra for CS replicate 3 
#' 
#' @param outputPath output folder path where all the text files will be writtem
#'
#' @return Six peaklist files each belonging to one replicate
#' 
#' @note: This script requires packages maldiquant, MALDIquantForeign and the custom R script exportMassSpectrumObjects.R
#'
#' @author Purva Kulkarni
#' 

pairwiseAnalysis-CJCS <- function(CJ_Replicate_1_path,CJ_Replicate_2_path,CJ_Replicate_3_path,CS_Replicate_1_path,CS_Replicate_2_path, CS_Replicate_3_path)
{
  # Generate list of mass spectrum objects for replicates of CJ
  CJ2_1b <- exportMassSpectrumObjects(CJ_Replicate_1_path)
  CJ2_2b <- exportMassSpectrumObjects(CJ_Replicate_2_path)
  CJ2_3 <- exportMassSpectrumObjects(CJ_Replicate_3_path)
  
  # Generate list of mass spectrum objects for replicates of CJ
  CS2_1 <- massFeatures(CS_Replicate_1_path)
  CS2_2 <- massFeatures(CS_Replicate_2_path)
  CS2_3 <- massFeatures(CS_Replicate_3_path)
  
  ## Average each mass spectrum list to create a single mass spectrum object (1 mass spectrum object per list)
  avgCJ2_1b <- averageMassSpectra(CJ2_1b, method = "mean")
  avgCJ2_2b <- averageMassSpectra(CJ2_2b, method = "mean")
  avgCJ2_3 <- averageMassSpectra(CJ2_3, method = "mean")
  avgCS2_1 <- averageMassSpectra(CS2_1, method = "mean")
  avgCS2_2 <- averageMassSpectra(CS2_2, method = "mean")
  avgCS2_3 <- averageMassSpectra(CS2_3, method = "mean")
  
  ## Create a list of the averaged mass spectrum objects
  LAESIspectra <- c(avgCJ2_1b, avgCJ2_2b, avgCJ2_3, avgCS2_1, avgCS2_2, avgCS2_3)
  
  ## Perform data preprocessing on this list to generate a list of filtered mass peak objects
  LAESIspectraSmoothed <- smoothIntensity(LAESIspectra, method="SavitzkyGolay", halfWindowSize=10)
  LAESIspectraSmoothedBaselineCorrected <- removeBaseline(LAESIspectraSmoothed, method="SNIP",iterations=100)
  LAESIspectraSmoothedBaselineCorrectedNormalized <- calibrateIntensity(LAESIspectraSmoothedBaselineCorrected, method="TIC")
  LAESIspectraSmoothedBaselineCorrectedNormalizedAligned <- alignSpectra(LAESIspectraSmoothedBaselineCorrectedNormalized)
  LAESIPeaks <- detectPeaks(LAESIspectraSmoothedBaselineCorrectedNormalizedAligned, SNR=4, halfWindowSize=10,  method="SuperSmoother")
  LAESIPeaksBinned <- binPeaks(LAESIPeaks, tolerance = 0.02)
  LAESIPeaksBinnedFILTERED <- filterPeaks(LAESIPeaksBinned, minFrequency=0.25)
  
  ## Export each peaklist in a separate file
  exportTab(LAESIPeaksBinnedFILTERED[[1]], file="LAESIPeaksBinnedFILTEREDCJ2_1b.tab")
  exportTab(LAESIPeaksBinnedFILTERED[[2]], file="LAESIPeaksBinnedFILTEREDCJ2_2b.tab")
  exportTab(LAESIPeaksBinnedFILTERED[[3]], file="LAESIPeaksBinnedFILTEREDCJ2_3.tab")
  exportTab(LAESIPeaksBinnedFILTERED[[4]], file="LAESIPeaksBinnedFILTEREDCS2_1.tab")
  exportTab(LAESIPeaksBinnedFILTERED[[5]], file="LAESIPeaksBinnedFILTEREDCS2_2.tab")
  exportTab(LAESIPeaksBinnedFILTERED[[6]], file="LAESIPeaksBinnedFILTEREDCS2_3.tab")
}

