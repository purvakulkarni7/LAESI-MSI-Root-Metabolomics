#' generate_volcano_plot.R
#' This function takes a csv file as input that contains fold change values for a sample pair and generates a volcano plot using it. Apart from this, it also prints the number of significant metabolites for each subset.
#'                                                      
#' @param filePath path to the *.csv file that contains the fold change values. This file can be exported from metaboAnalyst (http://metaboanalyst.ca) ver 4.0 after a volcano plot has been generated. The file mainly 5 columns and looks like the following:
#' 	FC	log2(FC)	raw.pval	#NAME?
#'892.2366443	7.0393	2.8154	0.00037312	3.4282
#'837.1989314	8.3037	3.0537	0.00062777	3.2022
#'
#' @export Exports *.png file with volcano plot. Along with this the function would also print the number of significant metabolites for each subset. 
#' 
#' @note The files used as an input for this code are:
#' Volcanoplot_Metaboanalyst_CSCJ.csv
#' Volcanoplot_Metaboanalyst_GMGP.csv
#' These input files are exported from metaboanalyst. Code snippets from metaboanalyst are saved in the files:
#' Metaboanalyst_code_pairwise_analysis_CSCJ.R
#' Metaboanalyst_code_pairwise_analysis_GMGP.R
#'
#' @author Purva Kulkarni

generate_volcano_plot <- function(filePath)
{
  
  # Read the metaboanalyst *.csv file 
  foldChangeTable <- read.table(filePath, header = TRUE, sep = ",") 
  
  # Make a basic volcano plot and export it as a *.png file
  par(bg=NA)
  with(foldChangeTable, plot(log2.FC., -log10(raw.pval), pch=20, col = "gray43", main="", xlab = "log2(FC)", ylab ="-log10(p)"))
  
  # Add lines
  abline(h = 1.0, col = "deepskyblue", lty = 2, lwd = 1)
  abline(v = c(-1,1), col = "deepskyblue", lty = 2, lwd = 1)
  
  # Add colored points
  with(subset(foldChangeTable, log2.FC.< -1 & X.log10.p.>1.0), points(log2.FC., -log10(raw.pval), pch=20, col="red"))
  with(subset(foldChangeTable, log2.FC.> 1 & X.log10.p.>1.0), points(log2.FC., -log10(raw.pval), pch=20, col="green3"))
  dev.copy(png,'volcano_plot.png')
  dev.off()
  
  # Print significant metabolites in each subset
  print(dim(foldChangeTable))
  print(dim(subset(foldChangeTable, log2.FC.< -1 & X.log10.p.>1.0)))
  print(dim(subset(foldChangeTable, log2.FC.> 1 & X.log10.p.>1.0)))
}