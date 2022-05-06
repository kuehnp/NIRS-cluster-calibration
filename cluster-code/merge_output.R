suppressPackageStartupMessages(library(foreach)) 						#suppress library notifications in the log

args <- commandArgs(trailingOnly=TRUE) 									#receive argument from finished 1st stage of array_job_wrapper.sh
input_dir <- args[1]													#use the first argument as input directory
input_file=list.files(input_dir,pattern="\\.rds$",full.names=TRUE)		#define input_file directory
resno<-args[2] 															#number of bootstrap results to bind together

finalMatrix <- foreach(i=1:resno, .combine=rbind) %do%{ # write the finalmatrix object by passing each optimizePLS.rds file to foreach, let it combine them sequentially(i.e. %do%)
	Optimized_Par<-readRDS(file=input_file)
	return(Optimized_Par) 
}

saveRDS(finalMatrix, file=file.path(input_dir,"finalmatrix.rds"))