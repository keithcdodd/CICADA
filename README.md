# CICADA
Comprehensive Independent Component Analysis Denoising Assistant (for resting and task-based functional magnetic resonance analysis)

CICADA can be broken down into three parts:
(1) Automatic CICADA: Automatically denoising fMRI images
(2) Manual CICADA: Save time with manual ICA denoising by only having to examine and adjust a small subset of ICs through CICADA
(3) Group CICADA: Run group-level quality control analyses and processing

Some details on installation and use can be seen below. More details can be found in the word document "CICADA_Details". Examples and more details on executing CICADA can also be found in the "example_CICADA_flow" directory.

# Installation
CICADA requires FSL (tested on 6.0) and Matlab (tested on 2022a). They must be configured to work together (e.g., call_fsl() in matlab calls fsl functions appropriately). Otherwise, this can be just be downloaded as a folder for the user, and the whole folder should be added to Matlab path (e.g., addpath(genpath()) to include all folders). See included word document for more installation and testing details. You can also find more information on workflow and testing the functions in the "example_CICADA_flow" directory.

# Requirements/Inputs
CICADA can work with any preprocessed data (warped to adult MNI space, but ideally not yet smoothed) but is well optimized for running on fMRIPrep outputs (tested on outputs from 20.2.6). CICADA can take advantage of the subject specific segmented probability masks (e.g., from fMRIPrep/FreeSurfer) and uses confound outputs (e.g., from fMRIPrep) to aid in selection of noise components. If fMRIPrep is not run, CICADA can be run without it (by calling the Auto_CICADA and Manual_CICADA functions). In its simplest form, one can use the template space's GM, WM, and CSF masks to aid in noise identification selection. CICADA was designed to be robust and effective for both resting-state and task-based fMRI datasets. It also works well with data with both low or high motion. See included word document ("CICADA_Details") and example(s) ("example_CICADA_flow") directory for more details.

# General Methods
Automatic CICADA:
First, CICADA creates regional masks of the functional data, and generates ICs with FSL's MELODIC. Next, relevant calculations are made from the regional masks, timeseries, and powerspectrum. These calculations are used to resort the ICs based on relative neural signal probability. Then, CICADA loops through a subset of these resorted ICs and uses the previous calculations to determine if the ICs should be labeled as signal or noise. Tagging of ICs of various signal-like and noise-like classifications is based on 3-group k-means calculations. CICADA then performs nonaggressive denoising of the noise ICs with FSL's fsl_regfilt, and then performs QC analyses. QC analyses include comparing CICADA denoising to standard 8 parameter denoising (6 motion parameters, white matter, and CSF) in a number of plots all representing common noise sources in fMRI. See included word document for more details.

Manual CICADA:
Following user-adjusted IC labeling, Manual CICADA reruns fsl_regfilt and QC analyses.

Group CICADA:
Runs similar QC analyses as Automatic CICADA, but on the group level. Also helps identify data outliers and prepares the data for easy analysis.

# Example Walkthrough
Automatic CICADA:
All relevant functions in CICADA can be called simply by using the included Auto_CICADA.m function. If the data was processed in fMRIPrep, it may be possible/easier to instead use the fmriprep_auto_CICADA.m function (which will, in turn, call the Auto_CICADA.m function). See included word document for more details. Examples of calling these functions can also be found in the example_CICADA_flow folder.

Manual CICADA:
First, adjust the Signal Label Column in the "IC_auto_checker.csv" (1 for true signal, 0 for true noise). Then re-save the file as "IC_manual_checker.csv". Now, run with Manual_CICADA.m or fmriprep_manual_CICADA.m. See included word document for more details.

Group CICADA:
Run cicada_group_qc function. See included word document for more details.

# Other Notes:
The MNI template that is used in fMRIPrep is provided. This can be utilized for more general analysis if you have warped to the same space, but do not have or want to use the subject specific segmentation.

# Citations
Publication is currently under review: Dodd K, McHugo M, Sarabia L, Wylie KP, Legget KT, Cornier MA, Tregellas JR. CICADA: An automated and flexible tool for comprehensive fMRI noise reduction. Imaging Neuroscience, 2024.

