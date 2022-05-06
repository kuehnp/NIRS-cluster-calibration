##### ----------- Read this section before you start --------------------- #####
#                                                                              #
# Script Version 6.3                                                           #
#                                                                              #
# This script covers the complete work flow                                    #
# from raw spectral files to predicted values in 12 steps:                     #
#                                                                              #
# Step 1  Define operation folders                                             #
# Step 2  Splicing and converting .asd files to .txt files                     #
# Step 3  Aggregating multiscans to single .txt files                          #
# Step 4  Import spectral data and set calibration parameters                  #
# Step 5  Creating a testset                                                   #
# Step 6  Export data to the cluster                                           #
#                                                                              #
# This script accesses the defined folders and writes on the disk.             #
# Making a backup of your all files before starting is highly recommended.     #
#                                                                              #
# Variables that require user input are in capital letters: VARIABLE ##        #
#                                                                              #
# A performance optimization and troubleshooting section                       #
# is included at the end of the file                                           #
#                                                                              #
################################################################################


rm(list=ls())     
gc()

# Loading Libraries
library(devtools)
#install_github(repo = "griffithdan/plantspec", force =TRUE) # manually install from GitHub (Requires Rtools and devtools)
#install_github(repo = "griffithdan/plantspecDB", force =TRUE) # manually install from GitHub (Requires Rtools and devtools)

library(plantspec)
library(plantspecDB)
library(Hmisc)
library(foreach)
library(doParallel)
library(htmlwidgets)
library(asdreader)
library(spectacles)

library(dplyr)
library(tidyr)
library(magrittr)
##OR##
library(tidyverse)


options(max.print=2000)


##### Step 1:  Define operation folders #####

SPECTRA_RAW          <- ("C:/NIRS/Spectra_for_Calibration_Raw")     ## Path that contains raw unprocessed spectral .asd files ()
SPECTRA_SPLICED      <- ("C:/NIRS/Spectra_for_Calibration_Spliced") ## spliced spectra will be written to this path
SPECTRA_CALIBRATION  <- ("C:/NIRS/Spectra_for_Calibration")         ## averaged spectra will be written to this path

SPECTRA_PREDICT_RAW  <- ("C:/NIRS/Spectra_to_be_Predicted_Raw")     ## Path that contains spectral files that 
SPECTRA_PREDICT      <- ("C:/NIRS/Spectra_to_be_Predicted")         ## Path that contains spectral files that 

MODELS_PATH          <- ("C:/NIRS/Prediction_Models")               ## Prediction models will be written to this path
RESULTS_PATH         <- ("C:/NIRS/Results")                         ## Predicted values will be written to this path



##### Step 2:  Splicing and converting asd files to txt files #####
# only necessary if you have ASD FieldSpec4 data files, otherwise skip this step
# requirement: One folder that contains all raw spectral files to convert "/SPECTRA_RAW/"

SPLICING = TRUE   ## Spicing is recommended for instruments with more than one sensor like ASD FieldSpec4
PLOTTING = TRUE   ## Plotting is recommended to visually quick-check the splicing results
FILES_LIST <- list.files(path = SPECTRA_RAW)

run <- length(FILES_LIST)

# Start loop
for (i in 1:run) {
  
  # Read spectrum
  MySpectrum <- get_spectra(paste0(SPECTRA_RAW, "/", FILES_LIST[i]))
  
  if (SPLICING == FALSE) {
    # Extracting spectral information
    Spectral_Numbers <- dimnames(MySpectrum)[[2]]
    Spectral_Data    <- as.numeric(MySpectrum)
    # Creating dataframe for output  
    OUTPUT_FRAME <- data.frame(Spectral_Numbers, Spectral_Data)
    FILENAME <- paste0(FILES_LIST[i], ".txt")
    
    # Plot Spectrum  
    if (PLOTTING ==TRUE) {plot(OUTPUT_FRAME, main=FILENAME)}
    
    # Export spectral file    
    write.table(OUTPUT_FRAME, paste0(SPECTRA_SPLICED, "/", FILENAME) , row.names = FALSE, col.names=FALSE)
  }
  
  if (SPLICING == TRUE) {
    # Extracting spectral information
    my.wl <- as.numeric(dimnames(MySpectrum)[[2]])
    my.id <- FILES_LIST[i]
    my.nir <- as.numeric(MySpectrum)
    
    # Creating a Spectra Object for splicing
    my.s <- Spectra(wl = my.wl, nir = my.nir, id = my.id)
    
    # Splicing
    my.s_spliced <- splice(my.s, locations = list(c(750, 1000), c(1950, 1800))) ## splicing points are configured for ASD FieldSpec4
    
    # Creating dataframe for output
    Spectral_Numbers <- my.s_spliced@wl
    Spectral_Data    <- as.numeric(my.s_spliced@nir[1,])
    OUTPUT_FRAME <- data.frame(Spectral_Numbers, Spectral_Data)
    FILENAME <- paste0(FILES_LIST[i], ".txt")
    
    # Plot Spectrum # splicing points are marked in blue 
    if (PLOTTING ==TRUE) {plot(OUTPUT_FRAME, main=FILENAME, axes=FALSE)
      axis(side=1, at=seq(350, 2500, by=100))
      axis(side=2, at=seq(0, 1, by=0.1))
      box()
      abline(v=1000, col="blue")
      abline(v=1800, col="blue")
    }
    
    # Export spectral file    
    write.table(OUTPUT_FRAME, paste0(SPECTRA_SPLICED, "/", FILENAME) , row.names = FALSE, col.names=FALSE)
  }
}

# Check for completion
OUTPUT_LIST <- list.files(path = SPECTRA_SPLICED)
OUTPUT_LIST <- gsub(".txt", "", OUTPUT_LIST)
identical(FILES_LIST, OUTPUT_LIST)  # will give TRUE if no error occurred 

# Generate a list of missing files
#MATCH_LIST <- match(FILES_LIST, OUTPUT_LIST)
#LOST_INDEX <- which(is.na(MATCH_LIST))
#MISSING_FILES <- FILES_LIST[LOST_INDEX]
#MISSING_FILES
#write.csv(MISSING_FILES, file="missing.csv")


##### Step 3:  Aggregating multiscans to single files #####
# only required if you have more than one spectral file per sample, otherwise skip this step
# requirement: - One folder that contains all spectral files to convert "/SPECTRA_SPLICED/"
#              - One csv file that contains one entry for each spectral file an information to which sample the file belongs   

# Read spectra
SPECTRA <- readSpectra(filelist = list.files(SPECTRA_SPLICED, full.names = T), wave_unit = "nm", measurement_unit = "reflectance", sep="", header=F)
#plot.spectra.matrix(SPECTRA)

# Read reference data
# The reference file must include a File and Sample variable
SAMPLE_REFERENCE <- read.csv("C:/NIRS/samples.csv", header = T, sep=";")
head(SAMPLE_REFERENCE)
length(SAMPLE_REFERENCE$Filename)
table(table(SAMPLE_REFERENCE$Filename))

# Check if the Reference matches the Spectra    
identical(SAMPLE_REFERENCE$Filename, dimnames(SPECTRA)[[1]])
length(SAMPLE_REFERENCE$Filename)
length(dimnames(SPECTRA)[[1]])

#SAMPLE_REFERENCE needs reordering
SAMPLE_REFERENCE<-SAMPLE_REFERENCE[match(dimnames(SPECTRA)[[1]], SAMPLE_REFERENCE$Filename),]

table(dimnames(SPECTRA)[[1]] %in% SAMPLE_REFERENCE$Filename )


testcheck<-as.data.frame(cbind( dimnames(SPECTRA)[[1]],SAMPLE_REFERENCE$Filename))
table(testcheck)
identical(dimnames(SPECTRA)[[1]], SAMPLE_REFERENCE$Filename)


# Averaging spectra on $Sample level  
SPECTRA <- as.spectra.matrix(SPECTRA)
SPECTRA_AV  <- averageSpectra(spec=SPECTRA, by=SAMPLE_REFERENCE$Sample)

# Export spectra    
writeSpectra(SPECTRA_AV, path=SPECTRA_CALIBRATION)


##### Step 4:  Import spectral data and set calibration parameters #####
# requirement: - One folder that contains all spectral files (1 spectral file = 1 sample, spliced and aggregated if necessary) "/SPECTRA_CALIBRATION/"
#              - One csv file that contains one entry for each spectral file with the according reference value
#              - Samples and reference values must be in the same order 

#Import Spectral Data
SPECTRA <- readSpectra(filelist = list.files(SPECTRA_CALIBRATION, full.names = T), wave_unit = "nm", measurement_unit = "reflectance")

# SPECTRA <- as.spectra.matrix(SPECTRA)

any(is.na(SPECTRA)) #FALSE

# Import Reference Dataset
DATA_REFERENCE <- read.csv("C:/NIRS/greenhouse_dry/traits.csv", sep=";", header=T)
head(DATA_REFERENCE) # check which reference data is available

#compare reference trait data to available spectra files
f<-list.files(SPECTRA_CALIBRATION, full.names = T)
#324, Reference_Traits.csv has 326 entries
refspectra<-as.data.frame(cbind(f))
refspectra$f
refspectra$samplecode<-word(refspectra$f, 5,sep = fixed("/")) %>% word (1, sep=fixed("."))
table(DATA_REFERENCE$Sample %in% refspectra$samplecode)
DATA_REFERENCE<-DATA_REFERENCE[match(refspectra$samplecode, DATA_REFERENCE$Sample),]


checktest<-as.data.frame(cbind(DATA_REFERENCE$Sample, refspectra$samplecode))


refspecc<- DATA_REFERENCE[,3:5]
table(table(refspecc$ID.No))

refspectra$samplecode %in% refspecc$SampleSpectra
identical(refspectra$samplecode,DATA_REFERENCE$Sample)
refspectra_comp <- merge(refspectra,refspecc,by.x="samplecode",by.y="SampleSpectra",all.y=T)
table(is.na.data.frame(refspectra_comp))

# Create new Variable for clarity
REFERENCE <- log(DATA_REFERENCE[,9]) ## Select the reference variable here
REFERENCE                        # double check if the correct reference variable is selected

qqnorm(REFERENCE, pch = 1, frame = FALSE)  # Q-q plot to check if the variable is normally distributed
qqline(REFERENCE, col = "red", lwd = 2)
hist(REFERENCE)                            # Histogram to check the distribution of the variable



# Set Parameter for region list
segments_number   <- 8     ## Meaningful values are 5 - 11      #    8 is a good solution for most applications. Increasing this value above 8 is not recommended. Doing so will result in high system saturation and system instability.
minimum_segment   <- 72    ## Meaningful values are 42 - 250    #   72 is a good solution for most applications.
minimum_frequence <- 350   ## Minimum frequency of spectrometer #  350 for ASD FieldSpec4 #  3800 for Bruker MPA.
maximum_frequence <- 2500  ## Maximum frequency of spectrometer # 2500 for ASD FieldSpec4 # 12400 for Bruker MPA.
wavesequence <- seq(from =minimum_frequence+minimum_segment, to=maximum_frequence-minimum_segment, by=1)

##### Step 5:  Creating a test-set #####
# Select a training set for model fitting. Will give a TRUE/FALSE vector.
# specifying TRUE for training/calibration data and FALSE for test/validation set data.
# p defines the percentage of VALIDATION samples: 0.3 will give 70% calibration samples.

trainingset <- !(subdivideDataset(spectra = SPECTRA, 
                                  method = "PCAKS",             ## Validation Method       # KS or PCAKS are good solutions for most applications.
                                  component = REFERENCE, 
                                  p = 0.5,                       ## Size of Validation Set  # Meaningful values are 0.3 - 0.7 # 0.5 is a good solution for most applications.
                                  type = "validation",
                                  output="logical"))

save.image(file="C:/NIRS/prep_workspace.rda")
load(file="C:/NIRS/prep_workspace.rda")
