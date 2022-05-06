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
# Step 6  Run calibration                                                      #
# Step 7  Select model                                                         #
# Step 8  Run validation                                                       #
# Step 9  Inspect model parameters                                             #
# Step 10 Outlier detection                                                    #
# Step 11 Export Model                                                         #
# Step 12 Apply Model(s)                                                       #
#                                                                              #
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

##### Step 1:  Define operation folders #####

SPECTRA_RAW          <- ("C:/NIRS/greenhouse_dry/Spectra_for_Calibration_Raw")     ## Path that contains raw unprocessed spectral .asd files ()
SPECTRA_SPLICED      <- ("C:/NIRS/greenhouse_dry/Spectra_for_Calibration_Spliced") ## spliced spectra will be written to this path
SPECTRA_CALIBRATION  <- ("C:/NIRS/greenhouse_dry/Spectra_for_Calibration")         ## averaged spectra will be written to this path

SPECTRA_PREDICT_RAW  <- ("C:/NIRS/greenhouse_dry/Spectra_to_be_Predicted_Raw")     ## Path that contains spectral files that 
SPECTRA_PREDICT      <- ("C:/NIRS/greenhouse_dry/Spectra_to_be_Predicted")         ## Path that contains spectral files that 

MODELS_PATH          <- ("C:/NIRS/greenhouse_dry/Prediction_Models")               ## Prediction models will be written to this path
RESULTS_PATH         <- ("C:/NIRS/greenhouse_dry/Results")                         ## Predicted values will be written to this path



##### Step 4:  Import spectral data and set calibration parameters #####
# requirement: - One folder that contains all spectral files (1 spectral file = 1 sample, spliced and aggregated if necessary) "/SPECTRA_CALIBRATION/"
#              - One csv file that contains one entry for each spectral file with the according reference value
#              - Samples and reference values must be in the same order 

#Import Spectral Data
SPECTRA <- readSpectra(filelist = list.files(SPECTRA_CALIBRATION, full.names = T), wave_unit = "nm", measurement_unit = "reflectance")

finalMatrix<-readRDS("C:/NIRS/greenhouse_dry/finalmatrix.rds")
finalMatrix<-readRDS("C:/NIRS/greenhouse_dry/finalmatrix_ghd.rds")
finalMatrix<-readRDS("C:/NIRS/greenhouse_dry/finalmatrix_50.rds")

##### Step 7:  Select model #####
# Select best model (here: lowest RMSE)
run <- dim(finalMatrix)

#use this if you chose test-set validation
best_RMSE <- finalMatrix[1,1][[1]]$RMSEP[1]   
for (i in 1:run[1]){
  if (finalMatrix[i,1][[1]]$RMSEP[1] <= best_RMSE){
    best_RMSE <- finalMatrix[i,1][[1]]$RMSEP[1]
    
    a <- finalMatrix[i,1][[1]]
    b <- finalMatrix[i,2][[1]]
    c <- finalMatrix[i,3][[1]]
    Selected_Model <- list(a,b,c)
    names(Selected_Model) <- c("optimization_results", "param_subsets", "param_preproc")
  }
}

#use this if you chose cross validation
'best_RMSE <- finalMatrix[1,1][[1]]$RMSECV[1]  #use this if you chose cross validation
for (i in 1:run[1])
{
  if (finalMatrix[i,1][[1]]$RMSECV[1] <= best_RMSE){
    best_RMSE <- finalMatrix[i,1][[1]]$RMSECV[1]
    
    a <- finalMatrix[i,1][[1]]
    b <- finalMatrix[i,2][[1]]
    c <- finalMatrix[i,3][[1]]
    Selected_Model <- list(a,b,c)
    names(Selected_Model) <- c("optimization_results", "param_subsets", "param_preproc")
  }
}'
  
  
  
  
##### Step 8:  Run validation #####
# Use Best Parameter for Model Calibration 
 
model_comparison<-as.data.frame(matrix(nrow=25,ncol=4)) %>% 'colnames<-'(c("R2_Cal","R2_Val","RMSEP","rank"))
  for (i in 1:25){
    tempmodcomp<-calibrate(component = REFERENCE, 
                           spectra = SPECTRA, 
                           optimal_params = Selected_Model, 
                           optimal_model = i,                 # use the best model. however, you might want to manually check the top 10 models and decide to chose a different one
                           validation = "testset",           
                           training_set = trainingset,        
                           max_comps = 12)
    model_comparison[i,1]<-tempmodcomp$R2_Cal
    model_comparison[i,2]<-tempmodcomp$R2_Val
    model_comparison[i,3]<-tempmodcomp$RMSEP
    model_comparison[i,4]<-tempmodcomp$rank
}

model_comparison_4<-model_comparison  
model_comparison_25<-model_comparison  
model_comparison_50<-model_comparison

   Best_Model <- calibrate(component = REFERENCE, 
                          spectra = SPECTRA, 
                          optimal_params = Selected_Model, 
                          optimal_model = 4,                 # use the best model. however, you might want to manually check the top 10 models and decide to chose a different one
                          validation = "testset",           
                          training_set = trainingset,        
                          max_comps = 12)
  
##### Step 9:  Inspect model parameters #####
  # Inspect Model Parameters
  Best_Model$R2_Cal #[1] 0.6610569
  Best_Model$R2_Val #[1] 0.6128125
  Best_Model$RMSEP
  Best_Model$rank   #[1] 8           
  Best_Model$regions
  Best_Model$preproc
  Best_Model$data$component
  
# View Model
  plot.PLScalibration(Best_Model, ncomp = Best_Model$rank)    # check calibration accuracy select value from Best_Model$rank here, or modify carefully
  # plot.PLScalibration(Best_Model, plottype = "validation")    # check different model ranks   
  
  dry_PERDICTIONS <- predictPLS(object= Best_Model, newdata = SPECTRA)  # Combined calibration and validation plot
  plot (as.numeric(dry_PERDICTIONS) ~ REFERENCE)                        #
  abline(lm((as.numeric(dry_PERDICTIONS) ~ REFERENCE)))                 # 
  abline(a=0, b=1, col="blue")                                      # adding x=y line for comparison
  text(-0.5,0.5,"R?=0.6743458, RMSEP= 0.1875247")
  
  summary(lm(as.numeric(dry_PERDICTIONS) ~ REFERENCE))$r.squared        # combined R2 of calibration and validation. For quick assessment only. Typically R2_Cal and R2_Val are better quality indicators.
  coef(lm(as.numeric(dry_PERDICTIONS) ~ REFERENCE))[1]

  as.data.frame(cbind(as.numeric(dry_PERDICTIONS),REFERENCE)) %>% ggplot(data=.,aes(x=REFERENCE, y=dry_PERDICTIONS))+
    geom_point(size=2)+
    theme_cowplot()+
    theme(panel.grid= element_line(color = "lightgrey",size=0.5), axis.text=element_text(size=16,face="bold"),
          axis.title=element_text(size=16, face="bold"))+
    geom_smooth(method=lm,se=F,size=2)+
    geom_abline(size=2)+
    lims(x=c(-1,0.75),y=c(-1,0.75))+
    labs(x="log leafN measured",y="log leafN predicted")+
    annotate("text",x=-0.5,y=0.5,label="R²=0.6610569, RMSEP= 0.1881825",size=7)  
  
  
  ##### Step 10: Outlier detection #####
  # Removing removing outliers is not recommended without further indication and 
  # has to be done manually by removing them from the INPUT_PATH and the Reference.csv file.
  
  DIFF <- abs(as.numeric(PERDICTIONS) - REFERENCE)
  FValue_DIFF <- DIFF
  DIFF_2 <- DIFF^2 
  
  Limit999 <- qf(.999, df1=1, df2=length(DIFF))
  Limit99  <- qf(.99, df1=1, df2=length(DIFF))
  Limit95  <- qf(.95, df1=1, df2=length(DIFF)) 
  Limit90  <- qf(.90, df1=1, df2=length(DIFF))
  
  for(i in 1:length(DIFF)){
    Value1 <- (length(DIFF) -1) * DIFF[i]^2
    Value2 <- (sum(DIFF_2) - DIFF[i])
    
    FValue_DIFF[i] <- Value1 / Value2}
  
  plot(FValue_DIFF)
  
  abline(h=Limit999,col="red")   # alpha = 0.999
  abline(h=Limit99, col="blue")  # alpha = 0.99
  abline(h=Limit95, col="green") # alpha = 0.95
  abline(h=Limit90, col="black") # alpha = 0.90
  
  outlier_index_999 <- (which(FValue_DIFF>Limit999))
  outlier_index_99  <- (which(FValue_DIFF>Limit99))
  outlier_index_95  <- (which(FValue_DIFF>Limit95))
  outlier_index_90  <- (which(FValue_DIFF>Limit90))
  
  FValue_DIFF[outlier_index_999]
  FValue_DIFF[outlier_index_99]
  FValue_DIFF[outlier_index_95]
  FValue_DIFF[outlier_index_90]
  
  # Names of outlier Spectra
  dimnames(SPECTRA)[[1]][outlier_index_999]
  dimnames(SPECTRA)[[1]][outlier_index_99]
  dimnames(SPECTRA)[[1]][outlier_index_95]
  dimnames(SPECTRA)[[1]][outlier_index_90]
  
  plotSpectra(SPECTRA)
  ##### Step 11: Export Model #####
  
  
  save(Best_Model, file = paste0(MODELS_PATH, "/", "Model_gh_leafN_dry.RData"))
  Best_model<-load(file = paste0(MODELS_PATH, "/", "Model_gh_leafN_dry.RData"))
  
##### compare best models from fresh and dry calibrations#####
table(names(dry_PERDICTIONS) %in% names(fresh_PERDICTIONS))
drp<-as.data.frame(dry_PERDICTIONS)
frp<-as.data.frame(fresh_PERDICTIONS)
drp$sample<-row.names(drp)
drp$runningno<-parse_number(drp$sample)
frp$sample<-row.names(frp)
frp$runningno<-parse_number(frp$sample)
drp<-drp[order(drp$runningno),]
frp<-frp[order(frp$runningno),]


drp<-drp[drp$runningno %in% frp$runningno,]
frp<-frp[frp$runningno %in% drp$runningno,]
drp$runningno %in% frp$runningno


cosine(frp$fresh_PERDICTIONS,drp$dry_PERDICTIONS)   #0.9121907
cor(frp$fresh_PERDICTIONS,drp$dry_PERDICTIONS)      #0.8173676
cov(frp$fresh_PERDICTIONS,drp$dry_PERDICTIONS,)     #0.05243014
boxplot(frp$fresh_PERDICTIONS,drp$dry_PERDICTIONS)



#####
  
  
  ##### Step 12 Apply prediction models to spectral files #####
  # requirement: - One folder that contains all spectral files on which the prediction models will be applied "/SPECTRA_PREDICT/"
  #              - One folder that contains all prediction models  "/MODELS_PATH/"
  #              - One folder where the results will be written "/RESULTS_PATH/"
  # It is possible to apply only one or multiple prediction models at once  
  
  
  library(devtools)
  library(plantspec)
  library(plantspecDB)
  library(Hmisc)
  
  # read spectral files
  SPECTRA <- readSpectra(filelist = list.files(SPECTRA_PREDICT, full.names = T), wave_unit = "nm", measurement_unit = "reflectance", sep="")
  
  plot.spectra.matrix(SPECTRA)
  
  MODELS_LIST <- list.files(path = MODELS_PATH)
  FILES_LIST <- list.files(path = SPECTRA_PREDICT) 
  
  Predictions_Matrix <- matrix(nrow = length(FILES_LIST), ncol = length(MODELS_LIST))
  colnames(Predictions_Matrix) <- as.list(MODELS_LIST)
  rownames(Predictions_Matrix) <- as.list(FILES_LIST)
  Run <- length(MODELS_LIST)
  
  for (i in 1:Run){
    load(paste0(MODELS_PATH,"/", MODELS_LIST[i]))
    PERDICTIONS <- predictPLS(object= Best_Model, newdata = SPECTRA)
    Predictions_Matrix[,i] <- as.numeric(PERDICTIONS)
  }
  hist(REFERENCE, ylim=c(0,100))
  hist(Predictions_Matrix, ylim=c(0,100))
  #predictions are heavily skewed towards the middle of the value range, lows and highs are reduced
  
  
  write.csv(Predictions_Matrix, file=paste0(RESULTS_PATH,"/gh_leafN.csv"), row.names = TRUE, na = "NA")
  
  ##### End of scipt #####
  
  
  #####-------------------------   performance optimization   -------------------------#####
  #                                                                                   
  #   The underlying idea of this script is to run multiple parallel iterations
  #   with different randomly chosen variables to find the best calibration model.
  #   Aim of this section is to demonstrate which influence which variable
  #   has on the calibration process. The performance of the script is primary
  #   dependent of the memory size (RAM), number of CPU cores and performance
  #   of the CPU cores. The influence of hard disk speed is small,
  #   influence of memory speed/timing is neglectable. The runtime of the
  #   script further depends greatly on the number of samples, the complexity
  #   of the wave sequence and the number of ranks of the calibration model.
  #   The key to optimize the calibration process is to optimally fit all model
  #   parameters considering the individual dataset to fully utilise the computer.
  #   It is important to note that more time and calculation power does not
  #   necessarily increase the model quality. However it might help to find
  #   the optimal model. Allocating more resources (Time and Power) to calibration
  #   will converge the found model slowly to the optimal model.
  #                                                                                   
  #                                      ***                                          
  #                                                                                   
  #   Calibration type (Test-Set vs Cross) and training set                           
  #                                                                                   
  #      trainingset <- !(subdivideDataset(spectra = SPECTRA,                         
  #                                  method = "random",                               
  #                                  component = REFERENCE,                           
  #                                  p = 0.5,                                         
  #                                  type = "validation",                             
  #                                  output="logical"))                               
  #                                                                                   
  #   In this script, test-set validation is significantly faster in calibration       
  #   and validation than cross-validation (approximately factor 5-10, depending       
  #   on the dataset). Additionally, there is a lower risk of accidentally overfitting 
  #   the model via test-set validation. However, cross-validation is a legit method   
  #   and an excellent solution for small sample sizes (n<80). The "subdivideDataset"       
  #   command divides the dataset in training set (calibration samples) and test       
  #   set (validation samples). The calibration set should account for the complete
  #   variation in yor dataset, as well as for the variation that is expected in future
  #   predictions. There is no consensus on which is the best selection   
  #   method for a test set and what is the optimal ratio of calibration to validation 
  #   samples are. Kennard Stones algorithm and PCA are commonly used (here 'KS'       
  #   and 'PCAKS'). Using 'random' is also legit, in particular it would be a good     
  #   idea to use multiple 'random' sets once you found the optimal model parameters   
  #   (see below) in order to test model stability. However it is advised to use       
  #   the same calibration for all traits (if you analyze multiple traits) and         
  #   to check if both calibration and validation set contain all of your species      
  #   (if your dataset contains multiple species). For medium sized datasets (n=150)   
  #   it is advised to have more calibration samples (for example 70% calibration,     
  #   30% validation). For very large datasets (n=500) it is also okay to have a         
  #   balanced ratio (50/50), if you are worried about model stability.
  #   Having more validation then calibration samples is uncommon, but it might
  #   be a good idea if you would like to check if your model is overfitted.                                                             
  #   
  #                                                                                   
  #   If you wish to use cross validation (not recommended for large datasets):                                             
  #     - in step 5 you need to define a testset even if you wish to use cross-validation.
  #                (for example use 'random' 50/50)                                   
  #     - in step 7 replace all "$RMSEP" by "$RMSECV"                                   
  #     - in step 8 delete the line "training_set = trainingset,"                       
  #     - in step 8 change "validation = "testset"," to "validation = "LOO","           
  #                                                                                   
  #                                      ***                                          
  #                                                                                   
  #                           Wave sequence and complexity                            
  #                                                                                   
  #   segments_number <- 7                                                            
  #   minimum_segment <- 72                                                           
  #   minimum_frequence <- 3800                                                       
  #   maximum_frequence <- 12400                                                      
  #                                                                                   
  #   The 'segments_number' is the number of spectral regions that are analysed       
  #   in every run. This variable has the highest impact on the resource demand.      
  #   Increasing this value will exponentially increase the memory demand of the      
  #   analysis. Keep in mind that running the script on multiple cpu cores will
  #   thereby also multiply the memory demand. If the computer runs out of memory
  #   it will use a swap file on the harddisk which is orders of magnitude slower.
  #   If it runs out of swap space the process might crash.
  #   If a system has limited memory but a strong CPU,       
  #   it is advised to reduce this value and run multiple models instead. If you      
  #   wish to run very complex models (9+) it is advised to run it on only a few      
  #   or single CPU core. When in doubt, it is advised to carefully monitor the       
  #   memory utilisation over the complete calibration process (memory demand will    
  #   increase in the process). Reducing the wave sequence by increasing the          
  #   'minimum_frequence' or decreasing the 'maximum_frequence' will in most          
  #   cases not decrease resource demand but might in some cases affect model         
  #   quality. For example, if there are noisy regions at the edge of the spectrum,   
  #   manually exclude them could increase model quality. However it is not           
  #   advised to remove spectral regions without indication.                          
  #                                                                                   
  #                                      ***                                          
  #                                                                                   
  #                               Parallel analysis                                   
  #                                                                                   
  #   cores=detectCores(logical = FALSE)                                              
  #   cl <- makeCluster(cores)                                                        
  #   registerDoParallel(cl)                                                          
  #   on.exit(stopCluster(cl))                                                        
  #   finalMatrix <- foreach(i=1:300, .combine=rbind) %dopar% { . }                   
  #                         
  #   The core of this script is the repeated parallel analysis of different 
  #   randomly drawn configurations of wave sequences. By manipulating these factors, 
  #   you define how many of these configurations are tested. The more configurations 
  #   you test, the higher is the likelihood that you will come close to an "optimal"
  #   configuration. The number of cores can be automatically detected 
  #   "cores=detectCores(logical = FALSE)" or manually set "cores<-5". The number 
  #   in the loop should be devisable by the number of cores that you use. For example
  #   if you want to test 300 configurations on 5 cores, put "foreach(i=1:300."
  #   and the script will run 60 iterations of 5 parallel runs. Before you start
  #   your main analysis it would be a good idea to make a test run with one iteration
  #   "cl <- makeCluster(5) .foreach(i=1:5." and check the time via System.time(). It 
  #   is advised that you run on 50 percent of the cores or up to cores-1 (if your 
  #   system has 6 cores, use 3-5). Especially for many iterations (>100) this will
  #   be faster than running on all cores.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
  #
  #                                      ***
  #
  #                                Model complexity
  #
  #   Optimized_Par <-  optimizePLS(component = REFERENCE,
  #                              spectra = SPECTRA, 
  #                              max_comps=12,
  #                              training_set = trainingset,
  #                              region_list = Region_List_A,
  #                              parallel = FALSE)
  #
  #   The 'max_comps' variable defines the highest possible model rank that is tested
  #   in the calibration process. It is okay to slightly overestimate this value, as
  #   the calibration process will most likely result in a model that has a lower rank
  #   then the 'max_comps'. However, if the resulting model rank is equal to 'max_comps'
  #   it might be a good idea to rerun the analyses with a higher 'max_comps' value. In
  #   this case it is advised to manually check all model ranks which should reveal a
  #   peak (meaning that increasing the model rank to a certain point should increase 
  #   model quality but beyond that, model quality should decrease), otherwise the model
  #   might be overfitted. The rank of a model is typically between 6 and 12 (but might
  #   range between 5 and 15, depending on the trait). This variable has high impact on
  #   the resource demand. Decreasing this value will increase calculation speed, but 
  #   might negatively affect model quality.
  #
  #                                      ***
  #
  #####------------------------     trouble shooting    ------------------------#####
  #
  #   Reading Spectral Files fails
  #
  #   The script uses a proprietary data type call spectra.matrix
  #   in some rare cases, the spectral information is not recognized.
  #   Run SPECTRA <- as.spectra.matrix(SPECTRA) to fix this
  #
  #   Error in names(x) <- value : 'names' attribute [3] must be the same length as the vector [2]
  #   this error occurs if the wrong separator is used. Use sep="" to clarify.
  #
  #    
  #
  #                                      ***
  #
  #   Creating a Test-Set crash or stuck
  #
  #   This process sometimes crashes or gets stuck in an endless loop. For example if 
  #   you try to create two even numbered dataset using an uneven number of samples or 
  #   in similar circumstances. Try using method = "random", or change from "p = 0.5",
  #   to "p = 0.49" for troubleshooting.
  #
  #                                      ***
  #
  #   Calibration error
  #
  #   The calibration process runs as part as a dplyr parallel pipeline, which makes
  #   tracing errors difficult. In case of a crash it is advised to run the code
  #   without the %dopar% {} to track the problem (set "i<-1"). Error "" occurs when
  #   the number of spectral files does not match the number of reference values.
  #   Double check your reference file and spectral files. 
  #
  #   Error "`$<-.data.frame`(`*tmp*`, "spectra", value = c(0, 0, 0, 0, 0,  : replacement has x rows, data has y"
  #   (while x and y are two different numerical values)
  #
  #   This error occurs when the number of calibration spectra does not match the number of reference values.
  #   This occures for exampe if you load files with header while they dont contain one.
  #
  #   Error in `$<-.data.frame`(`*tmp*`, "F", value = numeric(0)) : replacement has x rows, data has y
  #   (while x = 0 and y is a numerical value)
  #
  #   This error occurs when the wave sequence does not match the range of the spectral files. 
  #   Check if you selected the correct minimum and maximum frequency.
  #   Also manually check the spectral files for their spectral range.
  #   Also check if your files contain header information and wheather this is considered while loading the file.
  #
  #                                      ***
  #   Validation error
  #   
  #   An error occurs when the wrong variable type is addressed in validation 
  #   (i.e. when test set variables are used for cross validation and vice versa). 
  #   Check the calibration type section above for correct variable names.	
  #
  #   Bad Calibration Results
  #
  #   Run the optimization without the "%dopar% {.}" to manually observe the calibration events. 
  #   A typical symptom of this is, that the vast majority of the calibration events gives 
  #   rank 1 or 2 as the best model rank, even if you selected a higher "max_rank". This means 
  #   that the algorithm cannot find a good correlation between the spectra and the reference files.
  #   Doublecheck the reference data, are these physiologically meaningful data or is there a high 
  #   amount of outliers? Is it possible that there was a line skip in the reference file? It might 
  #   be a good idea to test a sub-dataset (for example only the first 20 files + reference). If this
  #   works try to enlarge this dataset until you find the outliers. Finally keep in mind that not 
  #   everything can be calibrated, there might just not be a signal of what you are searching for 
  #   in the spectra. Therefore it might be a good idea to start your analysis with a trait that is
  #   known to be predictable, for example LDMC for FieldSpec or nitrogen for Bruker MPA - this way 
  #   you can make sure that the process itself works.
  #
  #   Bad Validation Results
  #
  #   Typically, the validation R2 is slightly worse than the calibration R2. A difference of 
  #   0.05 to 0.1 is normal, a 0.2 difference is concerning (for example R2_cal=0.85; R2_val=0.65).
  #   If there is a big discrepancy between these values, this is typically a sign that the model 
  #   is overfitted. In this case you should have a look in the 'finalMatrix' if there is a model 
  #   on position 2, 3, 4 or so that has a better R2_Cal to R2_Val ratio, this might a better choice
  #   even if the RMSE is worse. Furthermore you should analyse the test-set. Are all species in cal
  #   and val set? Are the values evenly distributed? Consider increasing the number of validation
  #   samples. Consider using a different algorithm for sample selection.
  #
  ##### End of file #####
  