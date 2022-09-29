#!/bin/sh
# CADICA - Complete AutoDenoise - ICA: Can follow fMRIPrep.
# This is a script you can run if you already ran the first section up to the Matlab portion.
# check if aggressive or nonaggressive denoising! I do aggressive if we are doing "for sure noise dominated ICAs"

# 102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187
#  These are the variables to change (note: you may also need to make edits to the Matlab File!)
currsubjids=(102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187)
home="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME"
preproc="${home}/Preproc_ICA_rest"
derivatives="${preproc}/derivatives"
fmriprepdata="${derivatives}/fmriprep"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA" # your folder containing the scripts and such
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds

# Takes into account multiple sessions
for l in "${currsubjids[@]}"
do
  echo "Running for Subject ${l}"
  for j in "${sessids[@]}"
  do
  echo "Running for session ${j}"
  # go to current fmriprep subject folder

  # go to derivatives folder & make a folder for CADICA and subject folder
  cd ${derivatives}

  CADICAfol="${derivatives}/CADICA"
  cd ${CADICAfol}
  CADICAsubfol="${CADICAfol}/sub-${l}" # main subject folder to work within
  cd ${CADICAsubfol} # go into subject folder
  sessiondir="${CADICAsubfol}/ses-${j}"
  cd ${sessiondir} # if this doesn't work, then previous script was not executed correctly
  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  funcfile="${sessiondir}/funcfile.nii.gz" # functional file
  melfol="${sessiondir}/melodic"

  # might need to add the path to Matlab ahead of time. Have not found an easy and consistent way to do this yet in bash.
  cd ${sessiondir}
  #sessiondir_ml=\'${sessiondir}\'
  /Applications/MATLAB_R2022a.app/bin/matlab -batch 'cadica_cleaning_simpler'
  # You can look at the report from melodic alongside the ICA decisions from the matlab to see if it is doing what you think.
  # We found it most robust to look at the masks and ICA files in fslviewer over only looking at the report.

  # then you run fsl_regfilt for nonaggressive denoising! Like so:
  cd ${sessiondir}
  fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/CalculatedNoiseICs.csv)" \
  -d ${melfol}/melodic_mix -m ${funcmask} -a -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_bold.nii.gz"
  ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_bold.nii.gz"

  # Smooth with 6mm -> 6 / 2.3548 =2.548 to get effective gaussian kernel sigma
  fslmaths ${ICADenoised} -kernel gauss 2.548 -fmean -mul ${funcmask} "s_${ICADenoised}"

  # High pass filter (detrend). sigma = 1 / (2*f*TR), need to readd temporal mean that filter removes
  fslmaths "s_${ICADenoised}" -Tmean temporalmean.nii.gz
  fslmaths "s_${ICADenoised}" -bptf 31.25 -1 -add temporalmean.nii.gz -mul ${funcmask} "hpf_s_${ICADenoised}"

  done
  echo
  echo
done
