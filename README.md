# CADICA
Comprehensive Automatic Denoising Independent Component Analysis for rs-fMRI. This will step-by-step perform an automated version of ICA denoising. This takes into account common strategies used in manual denoising, and thus will help clean up both all types of physiologic and non-physiologic noise (including, but not limited to, motion).

# Installation
CADICA requires FSL and matlab or octave (tested on matlab 2022a). Otherwise, this can be just be downloaded as a folder for the user. It is highly recommended that fMRIPrep has also been run on the data of interest ahead of time.

# Requirements/Inputs
CADICA assumes fMRIPrep has previously been run on the data and CADICA is written to run smoothly with fMRIPrep output. CADICA takes advantage of the subject specific segmented probability masks (from fMRIPrep/FreeSurfer) and confound outputs from fMRIPrep to aid in selection of noise components. If fMRIPrep is not run, CADICA can be easily modified to run without it. In it's simplest form, one can use the template space's GM, WM, and CSF masks and comment out the use of the confounds to aid in noise identification selection. CADICA should remain effective and robust for most data even with these changes.

# Example Walkthrough
Will fill in later

# Other Notes:
MNI template that is used in fMRIPrep is provided. This can be utilized for more general analysis if you have warped to the same space, but do not have or want to use the subject specific segmentation.

# Citations
Publication in the works
