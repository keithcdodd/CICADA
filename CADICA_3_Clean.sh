#!/bin/sh
# CADICA - Comprehensive AutoDenoise - ICA: Can follow fMRIPrep.
# "Comprehensive" in the sense it takes into account non-bold physiologic and scanner-related noise alongside motion
# This is made to be run easily following fMRIPrep.
# This is modeled similarly to ALT (automated labeling tool) project with significant differences.
# For example, this uniquely takes advantage of subject specific anatomy masks of WM, CSF, and Edge masks in addition to GM.
# Uses FSL and Matlab.

# Decisions on Labelling ICs as Signal vs Noise is Based on the Following:
# (1) The significant clusters of an IC should lie more in Grey Matter than WM, CSF, or outside the brain.
# (2) The power of the IC should lie more in BOLD-related physiologic frequencies (0.01-0.1 Hz) than below or above
# (3) The IC timeseries should not show highly intense or numerous spikes (conservative cut-off)
# (4) The IC timeseries should not strongly correlate to confounds (global signal, WM, CSF, and all 6 motion parameters) (conservative cut-off)

# For the majority of data, (1) and (2) are sufficient while (3) and (4) are added safeguards that may pick up on other obvious missed noise
# For higher selection (less noise selected as signal), (1) and (2) can be modified from requiring "more of the signal" to requiring "a majority of the signal."
# but this might be too aggressive for particularly noisy data sets. A medium selection option is also available which selects only if the mean of the proportion
# of GM overlap plus BOLD frequency power is over 50%. The default is the less selective version, which then works well with aggressive denoising following.

# This script should be run following labeling. This script applies the following:
# (1) non-aggressive denoising
# (2) 6mm smoothing
# (3) Detrending (High Pass Filtering)

# Alternatively, if denoising in a different program (e.g., CONN), this script is not necessary. Just denoise with GLM of CADICA_noise_selection.mat, high pass filter and detrend.
# If doing non-aggressive denoising, consider trying to be more selective in signal labeling with, for example with medium selectivity setting in the Matlab labeling code.
# Results have been compared to Manual Labeling, ICA-AROMA, and ALT. It appears robust.
# Certain factors may perform better if adjusted by the user (as discussed in notes above).
# Check for changing factors in the code to fit your data. Such as TR, filenames/locations, and scan length.

######################### Variables to Change For Your Data As Needed #############################################
# curr subjects can be written as just the numbers as output in fMRIPrep. E.g. sub-102 is simply 102
# e.g. (001 002 003)
currsubjids=(102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187)
# where your CADICA data is held following the first two scripts
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA" # your folder containing the CADICA scripts and such
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds
cut_off_period="100" # 100 second period. This equates to 1/100 - 0.01 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
selection="medsel" # lowsel, medsel, highsel. lowsel labels more noise as signal. highsel labels more signal as noise. Highsel may be too selective for some data.
filt="nonaggressive" # nonaggressive or aggressive denoising, performed after selection. We suggest lowsel with agg, or medsel with nonagg.
#################################################################################################################

# Takes into account multiple sessions
for l in "${currsubjids[@]}"
do
  echo "Running for Subject ${l}"
  for j in "${sessids[@]}"
  do
  echo "  Running for session ${j}"
  # go to current fmriprep subject folder
  currsubjfol="${CADICAfol}/sub-${l}"
  sessiondir="${currsubjfol}/ses-${j}"
  cd ${sessiondir}
  funcfile="${sessiondir}/funcfile.nii.gz"
  funcmask="${sessiondir}/funcmask.nii.gz"
  melfol="${sessiondir}/melodic"

  # then you run fsl_regfilt for nonaggressive or aggressive denoising! You can select your level of selection as well. This is all set at the beginning of the script.
  if [[ "${filt}" == "aggressive" ]]; then
    echo "    Performing Filtering with fsl_regfilt with ${selection} and ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_${selection}_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -a -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_${selection}_agg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_${selection}_agg_bold.nii.gz"
  else
    echo "    Performing Filtering with fsl_regfilt with ${selection} and ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_${selection}_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_${selection}_nonagg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_${selection}_nonagg_bold.nii.gz"
  fi

  # Smooth with 6mm -> 6 / 2.3548 =2.548 to get effective gaussian kernel sigma
  echo "    Performing Smoothing With Gaussian 6mm kernel, with sigma of 2.548"
  fslmaths ${ICADenoised} -kernel gauss 2.548 -fmean -mul ${funcmask} "s_${ICADenoised}"

  # High pass filter (also detrends). sigma = 1 / (2*TR*cut_off_freq), and also readd temporal mean that filter removes
  # E.G. 1/(2*2*0.01)=25 for sigma cut off. We have to be fancy in bash because its math is limited.

  sigma_cutoff_interim=$((${cut_off_period}00/(2*${TR})))
  sigma_cutoff=$(echo "${sigma_cutoff_interim:0:2}.${sigma_cutoff_interim: -2}")
  echo "    Performing High Pass Filtering with Sigma of ${sigma_cutoff} for Cut Off Period of ${cut_off_period} seconds"
  fslmaths "s_${ICADenoised}" -Tmean temporalmean.nii.gz
  fslmaths "s_${ICADenoised}" -bptf "${sigma_cutoff}" -1 -add temporalmean.nii.gz -mul ${funcmask} "hpf_s_${ICADenoised}"
  done
done
