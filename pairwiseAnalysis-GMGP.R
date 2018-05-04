#' pairwiseAnalysis-GMGP.R
#' This function contains mass spectrmetry data analysis steps for the pairwise comaprative analysis of plant species GM and GP. The function parameters used in this function as customized as per the data.
#' 
#' @param GM_Replicate_1_path folder path containing 50 text files with mass spectra for GM replicate 1
#' @param GM_Replicate_2_path folder path containing 50 text files with mass spectra for GM replicate 2
#' @param GM_Replicate_3_path folder path containing 50 text files with mass spectra for GM replicate 3
#' 
#' @param GP_Replicate_1_path folder path containing 50 text files with mass spectra for GP replicate 1
#' @param GP_Replicate_2_path folder path containing 50 text files with mass spectra for GP replicate 2
#' @param GP_Replicate_3_path folder path containing 50 text files with mass spectra for GP replicate 3 
#' 
#' @param outputPath output folder path where all the text files will be writtem
#'
#' @return Six peaklist files each belonging to one replicate
#' 
#' @note: This script requires packages maldiquant, MALDIquantForeign and the custom R script exportMassSpectrumObjects.R
#'
#' @author Purva Kulkarni
#' 

pairwiseAnalysis-GMGP <- function(GM_Replicate_1_path,GM_Replicate_2_path,GM_Replicate_3_path,GP_Replicate_1_path,GP_Replicate_2_path, GP_Replicate_3_path)
{
  # Generate list of mass spectrum objects for replicates of GM
  GM1 <- massFeatures(GM_Replicate_1_path)
  GM2 <- massFeatures(GM_Replicate_2_path)
  GM3 <- massFeatures(GM_Replicate_3_path)
  
  # Generate list of mass spectrum objects for replicates of GP
  GP1 <- massFeatures(GP_Replicate_1_path)
  GP2 <- massFeatures(GP_Replicate_2_path)
  GP3 <- massFeatures(GP_Replicate_3_path)
  
  ## Average each mass spectrum list to create a single mass spectrum object (1 mass spectrum object per list)
  avgGM1 <- averageMassSpectra(GM1, method = "mean")
  avgGM2 <- averageMassSpectra(GM2, method = "mean")
  avgGM3 <- averageMassSpectra(GM3, method = "mean")
  avgGP1 <- averageMassSpectra(GP1, method = "mean")
  avgGP2 <- averageMassSpectra(GP2, method = "mean")
  avgGP3 <- averageMassSpectra(GP3, method = "mean")
  
  ## Create a list of the averaged mass spectrum objects
  LAESIspectraGMGP <- c(avgGM1, avgGM2, avgGM3, avgGP1, avgGP2, avgGP3)
  
  ## Perform data preprocessing on this list to generate a list of filtered mass peak objects
  LAESIspectraSmoothedGMGP <- smoothIntensity(LAESIspectraGMGP, method="SavitzkyGolay", halfWindowSize=10)
  LAESIspectraSmoothedBaselineCorrectedGMGP <- removeBaseline(LAESIspectraSmoothedGMGP, method="TopHat",halfWindowSize=5)
  LAESIspectraSmoothedBaselineCorrectedNormalizedGMGP <- calibrateIntensity(LAESIspectraSmoothedBaselineCorrectedGMGP, method="TIC")
  LAESIspectraSmoothedBaselineCorrectedNormalizedAlignedGMGP <- alignSpectra(LAESIspectraSmoothedBaselineCorrectedNormalizedGMGP)
  LAESIPeaksGMGP <- detectPeaks(LAESIspectraSmoothedBaselineCorrectedNormalizedAlignedGMGP, SNR=0.2, halfWindowSize=1,  method="MAD")
  LAESIPeaksBinnedGMGP <- binPeaks(LAESIPeaksGMGP, tolerance = 0.02)
  LAESIPeaksBinnedFILTEREDGMGP <- filterPeaks(LAESIPeaksBinnedGMGP, minFrequency=0.25)
  
  ## Export each peaklist in a separate file
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[1]], file="LAESIPeaksBinnedFILTEREDGM1.tab")
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[2]], file="LAESIPeaksBinnedFILTEREDGM2.tab")
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[3]], file="LAESIPeaksBinnedFILTEREDGM3.tab")
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[4]], file="LAESIPeaksBinnedFILTEREDGP1.tab")
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[5]], file="LAESIPeaksBinnedFILTEREDGP2.tab")
  exportTab(LAESIPeaksBinnedFILTEREDGMGP[[6]], file="LAESIPeaksBinnedFILTEREDGP3.tab")
}