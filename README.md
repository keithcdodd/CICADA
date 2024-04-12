# CICADA
Comprehensive Independent Component Analysis Denoising Assistant (for resting and task-based functional magnetic resonance analysis)

# Installation
CICADA requires FSL (tested on 6.0) and Matlab (tested on 2022a). They must be configured to work together (e.g., call_fsl() in matlab calls fsl functions appropriately). Otherwise, this can be just be downloaded as a folder for the user.

# Requirements/Inputs
CICADA can work with any preprocessed (warped to adult MNI space, but ideally not yet smoothed) but is well optimized for running on fMRIPrep outputs. CICADA can take advantage of the subject specific segmented probability masks (e.g., from fMRIPrep/FreeSurfer) and uses confound outputs (e.g., from fMRIPrep) to aid in selection of noise components. If fMRIPrep is not run, CICADA can be modified to run without it. In its simplest form, one can use the template space's GM, WM, and CSF masks to aid in noise identification selection. CICADA was designed to be robust and effective for both resting-state and task-based fMRI datasets. It also works well with data with both low or high motion.

# General Methods
Will fill in later

# Example Walkthrough
Will fill in later

# Other Notes:
The MNI template that is used in fMRIPrep is provided. This can be utilized for more general analysis if you have warped to the same space, but do not have or want to use the subject specific segmentation.

# Citations
Publication in the works
