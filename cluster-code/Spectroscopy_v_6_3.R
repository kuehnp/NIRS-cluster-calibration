### Calibration script ###
load("/home/user/prep_workspace.rda")

# Loading Libraries
suppressPackageStartupMessages(library(plantspec))
suppressPackageStartupMessages(library(plantspecDB))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(htmlwidgets))
suppressPackageStartupMessages(library(asdreader))
suppressPackageStartupMessages(library(spectacles))
options(max.print=2000)
sessionInfo()

args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output <- args[2]
seed <- args[3]

set.seed(seed)

##### Step 1:  Define operation folders #####

SPECTRA_CALIBRATION  <- input_dir         ## averaged spectra will be written to this path

##### Step 4:  Import spectral data and set calibration parameters #####
# requirement: - One folder that contains all spectral files (1 spectral file = 1 sample, spliced and aggregated if necessary) "/SPECTRA_CALIBRATION/"
#              - One csv file that contains one entry for each spectral file with the according reference value
#              - Samples and reference values must be in the same order 

#Import Spectral Data
SPECTRA <- readSpectra(filelist = list.files(SPECTRA_CALIBRATION, full.names = T), wave_unit = "nm", measurement_unit = "reflectance")

# SPECTRA <- as.spectra.matrix(SPECTRA
# Set Parameter for region list
segments_number   <- 8     ## Meaningful values are 5 - 11      #    8 is a good solution for most applications. Increasing this value above 8 is not recommended. Doing so will result in high system saturation and system instability.
minimum_segment   <- 72    ## Meaningful values are 42 - 250    #   72 is a good solution for most applications.
minimum_frequence <- 350   ## Minimum frequency of spectrometer #  350 for ASD FieldSpec4 #  3800 for Bruker MPA.
maximum_frequence <- 2500  ## Maximum frequency of spectrometer # 2500 for ASD FieldSpec4 # 12400 for Bruker MPA.
wavesequence <- seq(from =minimum_frequence+minimum_segment, to=maximum_frequence-minimum_segment, by=1)

##### Step 6:  Run calibration #####
# Select number of parallel tested iteration (number of parallel cores)
# Creating Cluster for parallel computing

print(paste("Processing begins" , Sys.time()))


# Create random sequence of NIR regions within the preset parameter
repeat{
	customsequence <- as.vector(sample(wavesequence, segments_number-1, replace = FALSE, prob = NULL))
	customsequence <- sort(customsequence, decreasing = FALSE)
	if (min(abs(diff(customsequence))) > minimum_segment) {break}
}
print(paste("End of randomization procedure" , Sys.time()))

# Assemble region list
region_list <- vector(mode = "list", length = segments_number)
for (i in 1:segments_number){
	if      (i == 1)              {region_list[[i]] = c(minimum_frequence, customsequence[i])   }
	else if (i == segments_number){region_list[[i]] = c(customsequence[i-1], maximum_frequence) }
	else                          {region_list[[i]] = c(customsequence[i-1], customsequence[i]) }}
print(paste("End internal FOR loop" , Sys.time()))


boxplot(region_list, horizontal = T)  

# Optimize Model Parameters 
Optimized_Par <-  optimizePLS(component = REFERENCE, 
spectra = SPECTRA, 
max_comps=12,                   ## meaningful values are 5-15, 12 is a good solution for most applications
region_list = region_list,
training_set = trainingset,     
parallel = FALSE)               # using parallel = TRUE is currently not supported
print(paste("Optimize Par" , Sys.time()))


saveRDS(Optimized_Par, file=output)