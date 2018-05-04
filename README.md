# LAESI-MSI Data Analysis for Untargeted Root Metabolomics
This reporsitroy contains all the R scripts, related inputs files and output files for comparative untargeted root metabolomics for two pairs of native and range expanding plant species.

| Pair | Native |Range expanding |
| :--------: | :--------: |:--------: |
| 1   | *Centuria jacea* (CJ)  |*Centuria stoebe* (CS)   |
| 2   | *Geranium molle* (GM)   |*Geranium pyrenaicum* (GP)   |

Export the mass spectrum files for every sample in indiviual `*.txt` files from Proteoplot software and apply the following scripts in a step-wise manner.

**1. exportMassSpectra.R**  
This function reads extracted mass spectra (from Proteoplot software) present in a text file and writes a new text file with only mz and intensity values, by removing the additional header lines. This happens to all the files in the folder path provided.

**2. exportMassSpectrumObjects.R**  
This function reads all the mass spectrum text files present in a folder and generates a list of mass spectum objects (an object type of the R package [MALDIquant](https://cran.r-project.org/web/packages/MALDIquant/index.html)) from it. The length of this list is equal to the number of mass spectrum files present in the folder.  

**3. pairwiseAnalysis-CJCS.R**    
This function contains mass spectrmetry data analysis steps for the pairwise comaprative analysis of plant species CJ and CS. The function parameters used in this function are customized as per the data.

**4. pairwiseAnalysis-GMGP.R**    
This function contains mass spectrmetry data analysis steps for the pairwise comaprative analysis of plant species GM and GP. The function parameters used in this function are customized as per the data.

**5. Metaboanalyst_code_MVA_CJCS_GMGP.R**    
This script contains R code exported from [Metaboanalyst](http://www.metaboanalyst.ca/) with all the parameters used for the the MVA analysis using mass feature matrix generated for all the four samples.

**6. Metaboanalyst_code_pairwise_analysis_CSCJ.R**  
This script contains R code exported from Metaboanalyst with all the parameters used for the the MVA analysis using mass feature matrix generated for pairwsie analysis of samples CS and CJ.

**7. Metaboanalyst_code_pairwise_analysis_GMGP.R**  
This script contains R code exported from Metaboanalyst with all the parameters used for the the MVA analysis using mass feature matrix generated for pairwsie analysis of samples GM and GP.

**8. generate_volcano_plot.R**    
Steps 3 and 4 would generate two feature matrices that can be used as an initial input for Metaboanalyst to generate a table of fold change values for mass features. This function takes a csv file as input that contains fold change values for a sample pair and generates a volcano plot using it. Apart from this, it also prints the number of significant metabolites for each subset.

**9. generate_box_plot.R**  
This function takes a csv file as input that contains normalized values of the sample concentrations for the replicates of a pair of plant species  (ex. CS and CJ)
