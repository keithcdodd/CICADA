# CICADA
Comprehensive Independent Component Analysis Denoising Assistant (for resting and task-based functional magnetic resonance analysis)

# Installation
CICADA requires FSL and matlab or octave (tested on matlab 2022a). Otherwise, this can be just be downloaded as a folder for the user. It is highly recommended that fMRIPrep has also been run on the data of interest ahead of time.

# Requirements/Inputs
CICADA can work with any preprocessed (warped to adult MNI space, but ideally not yet smoothed) but is well optimized for running on fMRIPrep outputs. CICADA can take advantage of the subject specific segmented probability masks (e.g., from fMRIPrep/FreeSurfer) and uses confound outputs (e.g., from fMRIPrep) to aid in selection of noise components. If fMRIPrep is not run, CICADA can be easily modified to run without it. In its simplest form, one can use the template space's GM, WM, and CSF masks to aid in noise identification selection. CICADA should remain effective and robust for most data sets, including resting-state and task-based analyses.

# Example Walkthrough
Will fill in later

# Other Notes:
MNI template that is used in fMRIPrep is provided. This can be utilized for more general analysis if you have warped to the same space, but do not have or want to use the subject specific segmentation.

# Citations
Publication in the works
