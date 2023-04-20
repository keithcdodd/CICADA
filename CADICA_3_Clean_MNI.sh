#!/bin/sh
# This script should be run following labeling. This script applies the following:
# (1) Denoising (aggressive or nonaggresive based on labeling)
# (2) 6mm smoothing
# (3) Detrending (High Pass Filtering)

# Alternatively, if denoising in a different program (e.g., CONN), this script is not necessary. Just denoise with GLM of CADICA_noise_selection.mat, high pass filter and detrend.

# Check for changing factors in the code to fit your data. Such as TR, filenames/locations, and scan length.

######################### Variables to Change For Your Data As Needed #############################################
# curr subjects can be written as just the numbers as output in fMRIPrep. E.g. sub-102 is simply 102
# e.g. (001 002 003)
currsubjids=(102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187)
# where your CADICA data is held following the first two scripts
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA_MNI_github" # your folder containing the CADICA scripts and such
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds
cut_off_period="125" # 125 second period. This equates to 1/100 - 0.008 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
filt="nonaggressive" # nonaggressive or aggressive denoising, performed after selection. We suggest lowsel with aggressive, or medsel/highsel with nonaggressive.
#################################################################################################################
echo
echo
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
    echo "    Performing Filtering with fsl_regfilt with ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_dist_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -a -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz"
  else
    echo "    Performing Filtering with fsl_regfilt with ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_dist_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz"
  fi

  # Smooth with 6mm -> 6 / 2.3548 =2.548 to get effective gaussian kernel sigma
  echo "    Performing Smoothing With Gaussian 6mm kernel, with sigma of 2.548"
  fslmaths ${ICADenoised} -kernel gauss 2.548 -fmean -mul ${funcmask} "s_${ICADenoised}"

  # High pass filter (also detrends). sigma = 1 / (2*TR*cut_off_freq), and also re-add temporal mean that filter removes
  # E.G. 1/(2*2*0.008)=31.25 for sigma cut off. We have to be fancy in bash because its math is limited.

  sigma_cutoff_interim=$((${cut_off_period}00/(2*${TR})))
  sigma_cutoff=$(echo "${sigma_cutoff_interim:0:2}.${sigma_cutoff_interim: -2}")
  echo "    Performing High Pass Filtering with Sigma of ${sigma_cutoff} for Cut Off Period of ${cut_off_period} seconds"
  fslmaths "s_${ICADenoised}" -Tmean temporalmean.nii.gz
  fslmaths "s_${ICADenoised}" -bptf "${sigma_cutoff}" -1 -add temporalmean.nii.gz -mul ${funcmask} "hpf_s_${ICADenoised}"
  done
done
