# CADICA
Complete Automatic Denoising Independent Component Analysis for rs-fMRI

# Installation
CADICA requires FSL and matlab or octave (tested on matlab 2022a). Otherwise, this can be just be downloaded as a folder for the user.

# Requirements/Inputs
CADICA assumes fMRIPrep has previously been run on the data and CADICA is written to run smoothly with fMRIPrep output. CADICA takes advantage of the subject specific segmented probability masks (from fMRIPrep/FreeSurfer) and confound outputs from fMRIPrep to aid in selection of noise components. If fMRIPrep is not run, CADICA can be easily modified to run without it. In it's simplest form, one can use the template space's GM, WM, and CSF masks and comment out the use of the confounds to aid in noise identification selection. CADICA should remain effective and robust for most data even with these changes.

# Example Walkthrough
Will fill in later

# Citations
Publication in the works
