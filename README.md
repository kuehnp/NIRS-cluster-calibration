# NIRS-cluster-calibration
An R-workflow utilizing the plantspec package to perform parallel calibrations of NIRS models for leaf traits on a computing cluster

This workflow was developed for the EVE - High-Performance Computing Cluster of the Helmholtz-Zentrum für Umweltforschung GmbH - UFZ (https://wiki.ufz.de/eve/). As of now only .asd files are supported. A basic familiarity with server computing tools is required. 

Many thanks to Tobias Proß, MLU Halle, for writing the initial R code, and to Christian Krause, iDiv Halle-Jena-Leipzig, for his support with adapting the code to the HPC cluster and parallelizing the calibration. 

## Requirements:
R Version 4.0 or higher

A HPC cluster running SLURM

Access to the HPC cluster command line (e.g. via PuTTY)

Access to the HPC Cluster workspace (e.g. via Filezilla)

A toolchain built for R (e.g. via Easybuild). 

Several R packages, especially plantspec (https://github.com/griffithdan/plantspec), asdreader (https://cran.r-project.org/web/packages/asdreader/) and Hmisc (https://cran.r-project.org/web/packages/Hmisc/), but see the R scripts for a complete list. 

## Data preparation
In this step the dataset will be prepared for calibration on a local PC. Input files consist of one folder containing NIRS Spectra and two datatables, *samples.csv* and and *traits.csv*. The *prep_cluster.R* will take this data, process it and output a handy *prep_workspace.rda* file to upload to the server. Extremely close attention should be paid to the order of the sample IDs in the dataframe and the SPECTRA dataset to make sure that they match and correspond correctly. 

## Cluster preparation
Here we set up our cluster pipeline, consisting of an *array_job_wrapper.sh*, which starts a number of *submit.sh* scripts, with each one starting a job running the *SPECTROSCOPY_v_6_3.R* script. Once the submit scripts are finished, the *array_job_wrapper.sh* launches *submit_merge_output.sh* which starts a job running the *merge_output.R*. This script pulls the results from each single calibration instance and feeds them into the merger script, which finally writes the results as a handy *finalmatrix.rds* file to be downloaded back to the local PC. 

The core is the *SPECTROSCOPY_v_6_3.R* script. Each instance of this script runs a calibration with randomized spectral regions, with the region set resulting in the best model passed on.

## Cluster calculation
Send the job off to the cluster and monitor it. To initiate the process we put 

*bash -x array_job_wrapper.sh /home/user/Spectra_for_Calibration 25*

into the bash command line. We can check the status of our job by using the *squeue* command. There will be two jobs running: one that is split up in a number of jobs (25 in our example), and one that is waiting for the other to finish so it can merge the results (signified by the label (Dependency) in the *squeue* output. Note that the different sub jobs will be distributed across many different computation nodes by design. 

If a mistake is made, the job can be cancelled by inputting *scancel [JOBID]*. You can cancel one part of the array job by inputting the specific subjob ID (e.g. *123456_2*) or the entire set of jobs by using the wrapper ID *123456*. 

## Data analysis
Download the *finalmatrix.rds* to the local PC. *results_cluster.R* will run an analysis and present a table with the 10 best models created during the calibration process. 
