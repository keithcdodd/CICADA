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
# 102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187
currsubjids=(108)
# where your CADICA data is held following the first two scripts
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA_MNI_github" # your folder containing the CADICA scripts and such
sessids=(01) # session numbering 01 02
TR="2" # repetition time of scan in seconds
low_cutoff_period="125" # 125 second period. This equates to 1/100 - 0.008 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
filt="nonaggressive" # nonaggressive or aggressive denoising, performed after selection. We suggest lowsel with aggressive, or medsel/highsel with nonaggressive.
anatmaskdir=
#################################################################################################################
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
  if [ ! -d "cleaned" ]
  then
    mkdir "cleaned"
  fi
  cleanedfol="${sessiondir}/cleaned"
  cd ${cleanedfol}

  funcfile="${sessiondir}/funcfile.nii.gz"
  funcmask="${sessiondir}/funcmask.nii.gz"
  melfol="${sessiondir}/melodic"
  anatmaskdir="${sessiondir}/anatmasks"

  # Anatomical Mask Import from before
  fslmaths "${anatmaskdir}/Edgeprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleanedfol}/EdgeOnly_mask.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleanedfol}/WMOnly_mask.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleanedfol}/CSFOnly_mask.nii.gz"

  Edgemask="${cleanedfol}/EdgeOnly_mask.nii.gz"
  WMmask="${cleanedfol}/WMOnly_mask.nii.gz"
  CSFmask="${cleanedfol}/CSFOnly_mask.nii.gz"

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

  # calculate timeseries for CSF, WM, Edge, and just low GM in general.
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -uthr 0.05 -bin -mul ${funcmask} "${anatmaskdir}/notGM_prop_final.nii.gz"
  fslmeants -i ${ICADenoised} -o CSF_timeseries.txt -m ${CSFmask}
  fslmeants -i ${ICADenoised} -o WM_timeseries.txt -m ${WMmask}
  fslmeants -i ${ICADenoised} -o Edge_timeseries.txt -m ${Edgemask}
  fslmeants -i ${ICADenoised} -o notGM_timeseries.txt -m "${anatmaskdir}/notGM_prop_final.nii.gz"
  fslmeants -i ${ICADenoised} -o Global_timeseries.txt -m ${funcmask}

  # Now combine all files to make design matrix in text form
  paste Edge_timeseries.txt WM_timeseries.txt CSF_timeseries.txt > designregressors.txt

  # convert to .mat with Text2Vest
  Text2Vest designregressors.txt designregressors.mat

  # Regress out CSF, WM, and Edge, calculate temporal mean ahead of time, and then demean when running since you did not include an intercept in the design matrix
  # might need to calculate temporal mean first as this may remove the temporal mean, and so you will add it back in (like you did for filtering after this)
  fslmaths "${ICADenoised}" -Tmean temporalmean.nii.gz
  fsl_glm -i ${ICADenoised} -d "designregressors.mat" --demean --out_res="3r_${ICADenoised}"
  fslmaths "3r_${ICADenoised}" -add temporalmean.nii.gz "3r_${ICADenoised}"

  # High pass filter (also detrends). sigma = 1 / (2*TR*cut_off_freq), and also re-add temporal mean that filter removes
  # E.G. 1/(2*2*0.008)=31.25 for sigma cut off. We have to be fancy in bash because its math is limited.
  # With a TR greater than 1.5s, we should likely only high pass filter.
  # low pass filtering not always recommended, especially for TR > 1.5s. FSL also does a poor job of low pass filtering. Could use AFNI.
  high_sigma_cutoff_interim=$((${low_cutoff_period}00/(2*${TR})))
  high_sigma_cutoff=$(echo "${high_sigma_cutoff_interim:0:2}.${high_sigma_cutoff_interim: -2}")

  echo "    Performing High Pass Filtering with Sigmas of ${high_sigma_cutoff} for Cut Off Period of ${low_cutoff_period} seconds"
  fslmaths "3r_${ICADenoised}" -Tmean temporalmean.nii.gz
  fslmaths "3r_${ICADenoised}" -bptf "${high_sigma_cutoff}" -1 -add temporalmean.nii.gz -mul ${funcmask} "hpf_3r_${ICADenoised}"

  # Smooth with 6mm -> 6 / 2.3548 =2.548 to get effective gaussian kernel sigma. Recommended (e.g., in AFNI) to smooth after filtering.
  echo "    Performing Smoothing With Gaussian 6mm kernel, with sigma of 2.548"
  fslmaths "hpf_3r_${ICADenoised}" -kernel gauss 2.548 -fmean -mul ${funcmask} "s_hpf_3r_${ICADenoised}"

  done
done
echo

# Now we have fully preprocessed functionals (assuming this is resting state), and now can move on to analyses. For example, seed-based, or dual regression of group ICA...
