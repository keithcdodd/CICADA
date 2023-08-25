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
currsubjids=(103)
# where your CADICA data is held following the first two scripts
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated"
taskid="rest"
CADICAfuncs="/Users/keithdodd/GitHub/CADICA/Newer" # your folder containing the CADICA scripts and such
sessids=(01) # session numbering 01 02
TR="2" # repetition time of scan in seconds
lowfreq_cutoff="0.008" # in Hz
highfreq_cutoff="0.15" # in Hz
blur="6" # mm gaussian blur
filt="nonaggressive" # nonaggressive or aggressive denoising, performed after selection. We suggest nonaggressive.
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

  # funcfile="${sessiondir}/funcfile.nii.gz"
  orig_funcfile="${sessiondir}/funcfile.nii.gz"
  funcfile="${sessiondir}/s_bp_funcfile.nii.gz"
  funcmask="${sessiondir}/funcmask.nii.gz"
  melfol="${sessiondir}/melodic"
  anatmaskdir="${sessiondir}/anatmasks"
  anatprob="${anatmaskdir}/Anatprob_resam.nii.gz"

  fslmaths "${anatprob}" -thr 0.5 -bin "${anatmaskdir}/anatmask_final.nii.gz"

  anatmask="${anatmaskdir}/anatmask_final.nii.gz"

  # Anatomical Mask Import from before
  fslmaths "${anatmaskdir}/Edge_prop.nii.gz" -sub "${anatmaskdir}/GMprob_resam.nii.gz" -thr 0.95 -bin "${cleanedfol}/Edge_selective_mask.nii.gz"
  fslmaths "${anatmaskdir}/WMCSF_boundary_prop.nii.gz" -sub "${anatmaskdir}/GMprob_resam.nii.gz" -thr 0.95 -bin "${cleanedfol}/WMCSF_selective_mask.nii.gz"
  fslmaths "${funcmask}" -sub "${anatprob}" -sub "${anatmaskdir}/Edge_prop.nii.gz" -thr 0.95 -bin "${cleanedfol}/Outbrain_noEdge_selective_mask.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleanedfol}/WM_selective_mask.nii.gz"
  fslmaths "${anatmaskdir}/CSFprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleanedfol}/CSF_selective_mask.nii.gz"

  Edgeselectivemask="${cleanedfol}/Edge_selective_mask.nii.gz"
  WMCSFselectivemask="${cleanedfol}/WMCSF_selective_mask.nii.gz"
  Outbrainselectivemask="${cleanedfol}/Outbrain_noEdge_selective_mask.nii.gz"
  WMselectivemask="${cleanedfol}/WM_selective_mask.nii.gz"
  CSFselectivemask="${cleanedfol}/CSF_selective_mask.nii.gz"

  # Calculate generally selective masks too, in case you want them as regressors:
  fslmaths "${anatmaskdir}/Edge_prop.nii.gz" -thr 0.95 -bin "${cleanedfol}/Edge_selective_mask.nii.gz"

  fslmaths "${anatmaskdir}/WMCSF_boundary_prop.nii.gz" -thr 0.95 -bin "${cleanedfol}/WMCSF_selective_prop_final.nii.gz"

  Edgemask="${anatmaskdir}/Edge_prop_final.nii.gz"
  Subepemask="${anatmaskdir}/Subepe_prop_final.nii.gz"
  CSFmask="${anatmaskdir}/CSF_prop_final.nii.gz"
  OutbrainNoEdgemask="${anatmaskdir}/Outbrain_noEdge_prop_final.nii.gz"

  # then you run fsl_regfilt for nonaggressive or aggressive denoising! You can select your level of selection as well. This is all set at the beginning of the script.
  if [[ "${filt}" == "aggressive" ]]; then
    echo "    Performing Filtering with fsl_regfilt with ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_dist_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -a -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz"
    fslmaths "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz" -mul ${anatmask} "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz"
  else
    echo "    Performing Filtering with fsl_regfilt with ${filt} Denoising."
    fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/Noise_dist_ICs.csv)" \
    -d ${melfol}/melodic_mix -m ${funcmask} -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz"
    fslmaths "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz" -mul ${anatmask} "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz"
    ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz"
  fi

  # calculate timeseries for CSF, WM, Edge, and just low GM in general.
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -uthr 0.05 -bin -mul ${funcmask} "${anatmaskdir}/notGM_prop_final.nii.gz"
  fslmeants -i ${ICADenoised} -o CSF_CADICA_timeseries.txt -m ${CSFselectivemask}
  fslmeants -i ${ICADenoised} -o WM_CADICA_timeseries.txt -m ${WMselectivemask}
  fslmeants -i ${ICADenoised} -o Edge_CADICA_timeseries.txt -m ${Edgeselectivemask}
  fslmeants -i ${ICADenoised} -o notGM_CADICA_timeseries.txt -m "${anatmaskdir}/notGM_prop_final.nii.gz"
  fslmeants -i ${ICADenoised} -o Global_CADICA_timeseries.txt -m ${funcmask}
  fslmeants -i "${funcfile}" -o EdgeOnly_timeseries.txt -m "${Edgeselectivemask}"
  fslmeants -i "${funcfile}" -o SubepeOnly_timeseries.txt -m "${WMCSFselectivemask}"
  fslmeants -i "${funcfile}" -o CSFOnly_timeseries.txt -m "${CSFselectivemask}"
  fslmeants -i "${funcfile}" -o OutbrainOnly_timeseries.txt -m "${Outbrainselectivemask}"


  # Now combine ICA Denoised regressor files to make design matrix in text form
  paste Edge_timeseries.txt WM_timeseries.txt CSF_timeseries.txt > designregressorsICA.txt

  # Also combine the funcfile regressors to test
  paste EdgeOnly_timeseries.txt SubepeOnly_timeseries.txt CSFOnly_timeseries.txt OutbrainOnly_timeseries.txt > designregressors.txt

  # compare to standard 9 Parameters
  fslmaths "${funcfile}" -mul "${funcmask}" -Tmean temporalmean.nii.gz
  3dTproject -input "${orig_funcfile}" -ort "${sessiondir}/9p_regressors.txt" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "s_bp_9p_funcfile.nii.gz"
  fslmaths "s_bp_9p_funcfile.nii.gz" -add temporalmean.nii.gz -mul "${funcmask}" "s_bp_9p_funcfile.nii.gz"

  # can compare to just regression of the 4 noise masks:
  fslmaths "${funcfile}" -mul "${funcmask}" -Tmean temporalmean.nii.gz
  3dTproject -input "${funcfile}" -ort "designregressors.txt" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "s_bp_4r_funcfile.nii.gz"
  fslmaths "s_bp_4r_funcfile.nii.gz" -add temporalmean.nii.gz -mul "${funcmask}" "s_bp_4r_funcfile.nii.gz"

  :'
  # Bandpass and blur only (main, best, and final version)
  fslmaths "${ICADenoised}" -Tmean temporalmean.nii.gz
  3dTproject -input "${ICADenoised}" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "s_bp_${ICADenoised}"
  fslmaths "s_bp_${ICADenoised}" -add temporalmean.nii.gz -mul "${anatmask}" "s_bp_${ICADenoised}"


  # Can compare to also regressing out CSF, WM, and Edge, calculate temporal mean ahead of time to add back in, also bandpass and blur
  fslmaths "${ICADenoised}" -Tmean temporalmean.nii.gz
  3dTproject -input "${ICADenoised}" -ort designregressorsICA.txt -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "s_bp_3r_${ICADenoised}"
  fslmaths "s_bp_3r_${ICADenoised}" -add temporalmean.nii.gz -mul "${anatmask}" "s_bp_3r_${ICADenoised}"

  # compare to just bandpass and smooth original data
  fslmaths "${funcfile}" -mul "${funcmask}" -Tmean temporalmean.nii.gz
  3dTproject -input "${funcfile}" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "s_bp_funcfile.nii.gz"
  fslmaths "s_bp_funcfile.nii.gz" -add temporalmean.nii.gz -mul "${funcmask}" "s_bp_funcfile.nii.gz"
  '

  # Create a QC folder
  QCdir="${sessiondir}/QC"
  # rm -rf "${QCdir}"
  # mkdir "${QCdir}"
  # mkdir "${QCdir}/parctimeseries"

  done
done
echo

# Now we have fully preprocessed functionals (assuming this is resting state), and now can move on to analyses. For example, seed-based, or dual regression of group ICA...
