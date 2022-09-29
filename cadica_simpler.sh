#!/bin/sh
# CADICA - Complete AutoDenoise - ICA: Can follow fMRIPrep.
# "Complete" in the sense it takes into account non-bold physiologic and scanner-related noise alongside motion (no need to nuissance regress WM or CSF afterwards like ICA-AROMA)
# This is made to be run easily following fMRIPrep.
# This is modeled similar to ALT (automated labeling tool) project,
# This uniquely takes advantage of subject specific anatomy masks of WM, CSF, and Edge masks in addition to GM
# Uses FSL and Matlab. Uses subject specific masks of edge voxels, GM, WM, and CSF, alongside powerspectrum analysis from FSL melodic_IC
# Uses these to calculate what ICA aspects to keep with nonaggressive Denoising
# Decisions are based on three simple physiologic principles:
# (1) The majority of the component should not lie on or outside the edge of the brain (otherwise, it is likely motion dominant)
# (2) The component signal within the brain should proportionally lie within GM significantly more than random chance
# (3) The power of the component signal should lie more within BOLD-related physiologic frequencies (0.01 - 0.1 Hz) than below or above
# At the end, 6mm smoothing and detrending (High Pass Filtering) is applied which can complete the denoising procedure
# Results have been compared to ICA-AROMA, ALT, and Manual Labeling. It appears robust.
# Certain factors may perform better if adjusted by the user.
# You may need to edit the Matlab file (e.g. change TR & scan length as needed) and add to Matlab path
# This assumes fmriPrep output, so you don't need to give specific inputs this way

# done: 102, first session of 103... continue on with a few more before running all
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
  currsubjfol="${fmriprepdata}/sub-${l}"
  cd ${currsubjfol}

  # look for mask files, fmriprep is not consistent yet with this placement
  if [ -d "${currsubjfol}/anat" ]
  then
      anatfol="${currsubjfol}/anat"
      anatfilestart="sub-${l}"
      # keep track of if there is more than one anat file to worry about
      multanats="0"
  elif  [ -d "${currsubjfol}/ses-01/anat" ]
  then
      # because anat file will only be in session 01 folder
      anatfol="${currsubjfol}/ses-01/anat"
      anatfilestart="sub-${l}_ses-01"
      multanats="1"
  else
      echo "Cannot find relevant anatomy folder"
      exit 1
  fi

  # find functional folder too!
  funcfol="${currsubjfol}/ses-${j}/func"

  # go to derivatives folder & make a folder for CADICA and subject folder
  cd ${derivatives}

  # make CADICA (Complete Auto Denoise - ICA) folder if we have not already
  if [ ! -d "${derivatives}/CADICA" ]
  then
    mkdir CADICA
  fi
  CADICAfol="${derivatives}/CADICA"
  cd ${CADICAfol}

  # make subject folder if we have not already
  if [ ! -d "${CADICAfol}/sub-${l}" ]
  then
    mkdir "sub-${l}"
  fi
  CADICAsubfol="${CADICAfol}/sub-${l}" # main subject folder to work within
  cd ${CADICAsubfol} # go into subject folder

  # delete old session folder if it exists, and make new one
  if [ -d "${CADICAsubfol}/ses-${j}" ]
  then
    rm -rf "${CADICAsubfol}/ses-${j}"
  fi
  mkdir "ses-${j}" # make session folder to put all functional files into
  sessiondir="${CADICAsubfol}/ses-${j}"
  cd ${sessiondir}

  # make session anat mask folder if we have not already
  if [ ! -d "${CADICAsubfol}/ses-${j}/anatmasks" ]
  then
    mkdir "anatmasks" # make a mask folder to put all anatomical masks of interest into
  fi
  anatmaskdir="${sessiondir}/anatmasks"

  # copy over session functional files and confounds file
  cp "${funcfol}/sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz" ${sessiondir}/funcmask.nii.gz # func brain mask
  cp "${funcfol}/sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz" ${sessiondir}/funcfile.nii.gz # func file
  cp "${funcfol}/sub-${l}_ses-${j}_task-${taskid}_desc-confounds_timeseries.tsv" ${sessiondir}/confounds_timeseries

  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  funcfile="${sessiondir}/funcfile.nii.gz" # functional file

  # copy over anatomical stuff
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz" ${anatmaskdir}/anatmaskorig.nii.gz # 1 is GM, 2 is WM, 3 is CSF
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_dseg.nii.gz" ${anatmaskdir}/anatseg.nii.gz
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz" ${anatmaskdir}/GMprobseg.nii.gz
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz" ${anatmaskdir}/WMprobseg.nii.gz
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz" ${anatmaskdir}/CSFprobseg.nii.gz
  # also copy over the original anatomy file to main folder, just for easier reference later as desired (e.g., in CONN)
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz" ${sessiondir}/anatfile.nii.gz

  # relabel it all for easier coding
  anatseg="${anatmaskdir}/anatseg.nii.gz" # 1 is GM, 2 is WM, 3 is CSF
  GMprob="${anatmaskdir}/GMprobseg.nii.gz"
  WMprob="${anatmaskdir}/WMprobseg.nii.gz"
  CSFprob="${anatmaskdir}/CSFprobseg.nii.gz"
  anatmask="${anatmaskdir}/anatmaskorig.nii.gz" # mask of anat file

  # calculate in the brain and out of the brain based on the probability segments
  fslmaths "${GMprob}" -add "${WMprob}" -add "${CSFprob}" "${anatmaskdir}/inbrainprob.nii.gz"
  # for out of brain, we need to resample funcmask to same space as others for now
  flirt -ref ${anatmask} -in ${funcmask} -out "${anatmaskdir}/funcmask_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  fslmaths "${anatmaskdir}/funcmask_resam.nii.gz" -add "${anatmask}" -bin -sub "${anatmaskdir}/inbrainprob.nii.gz" -thr 0.01 "${anatmaskdir}/outbrainprob.nii.gz"

  # relabel These
  inbrainprob="${anatmaskdir}/inbrainprob.nii.gz"
  outbrainprob="${anatmaskdir}/outbrainprob.nii.gz"

  # resample anatseg, anatmask, anatprobs so it is all in the same space
  flirt -ref ${funcmask} -in ${anatseg} -out "${anatmaskdir}/anatseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm -interp "nearestneighbour"
  flirt -ref ${funcmask} -in ${anatmask} -out "${anatmaskdir}/anatmask_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm -interp "nearestneighbour"
  flirt -ref ${funcmask} -in ${GMprob} -out "${anatmaskdir}/GMprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${WMprob} -out "${anatmaskdir}/WMprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${CSFprob} -out "${anatmaskdir}/CSFprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${inbrainprob} -out "${anatmaskdir}/inbrainprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${outbrainprob} -out "${anatmaskdir}/outbrainprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  fslmaths "${anatmaskdir}/outbrainprobseg_resam.nii.gz" -mul "${funcmask}" "${anatmaskdir}/outbrainprobseg_resam_masked.nii.gz"

  # relabel these
  anatbrainmask="${anatmaskdir}/anatmask_resam.nii.gz" # from anat mask originally resampled
  GMprob_resam="${anatmaskdir}/GMprobseg_resam.nii.gz"
  WMprob_resam="${anatmaskdir}/WMprobseg_resam.nii.gz"
  CSFprob_resam="${anatmaskdir}/CSFprobseg_resam.nii.gz"
  inbrainprob_resam="${anatmaskdir}/inbrainprobseg_resam.nii.gz"
  outbrainprob_resam="${anatmaskdir}/outbrainprobseg_resam_masked.nii.gz"

  # let's finally get to melodic IC decomposition
  cd ${sessiondir}
  # make session folder if we have not already
  if [ ! -d "${sessiondir}/melodic" ]
  then
    mkdir melodic # make melodic folder to run melodic in
  fi
  melfol="${sessiondir}/melodic"

  # Run Melodic!
  melodic --in="${funcfile}" --outdir=${melfol} --mask=${funcmask} --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}

  cd ${sessiondir}
  # For each ICA: smooth & threshold standard (gives clusters), multiply by masks, and get useful output (e.g. numvoxels)

  # Full volume
  fslmaths ${melfol}/melodic_IC.nii.gz -fmean -thr 3 -mul "${funcmask}" "${sessiondir}/fullvolICA_adj.nii.gz"
  fslstats -t "${sessiondir}/fullvolICA_adj.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > fullICmean.txt
  awk '{print $2}' tmp.txt > fullICnumvoxels.txt

  # GM voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${GMprob_resam}" "${sessiondir}/GMICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/GMICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > GMICmean.txt
  awk '{print $2}' tmp.txt > GMICnumvoxels.txt

  # WM voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${WMprob_resam}" "${sessiondir}/WMICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/WMICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > WMICmean.txt
  awk '{print $2}' tmp.txt > WMICnumvoxels.txt

  # CSF voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${CSFprob_resam}" "${sessiondir}/CSFICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/CSFICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > CSFICmean.txt
  awk '{print $2}' tmp.txt > CSFICnumvoxels.txt

  # inbrain voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${inbrainprob_resam}" "${sessiondir}/inbrainICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/inbrainICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > inbrainICmean.txt
  awk '{print $2}' tmp.txt > inbrainICnumvoxels.txt

  # outbrain voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${outbrainprob_resam}" "${sessiondir}/outbrainICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/outbrainICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > outbrainICmean.txt
  awk '{print $2}' tmp.txt > outbrainICnumvoxels.txt

  rm tmp.txt

  # Then you run the matlab file
  # might need to add the path to Matlab ahead of time. Have not found an easy and consistent way to do this yet in bash.
  cd ${sessiondir}
  #sessiondir_ml=\'${sessiondir}\'
  /Applications/MATLAB_R2022a.app/bin/matlab -batch 'cadica_cleaning_simpler'
  # You can look at the report from melodic alongside the ICA decisions from the matlab to see if it is doing what you think.
  # We found it most robust to look at the masks and ICA files in fslviewer over only looking at the report.

  # then you run fsl_regfilt for nonaggressive denoising! Like so:
  cd ${sessiondir}
  fsl_regfilt -i ${funcfile} -f "$(cat ${sessiondir}/CalculatedNoiseICs.csv)" \
  -d ${melfol}/melodic_mix -m ${funcmask} -o "sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_bold.nii.gz"
  ICADenoised="sub-${l}_ses-${j}_task-${taskid}_space-MNI152NLin2009cAsym_desc-CADICA_bold.nii.gz"

  # Smooth with 6mm -> 6 / 2.3548 =2.548 to get effective gaussian kernel sigma
  fslmaths ${ICADenoised} -kernel gauss 2.548 -fmean -mul ${funcmask} "s_${ICADenoised}"

  # High pass filter (detrend). sigma = 1 / (2*f*TR), need to readd temporal mean that filter removes
  fslmaths "s_${ICADenoised}" -Tmean temporalmean.nii.gz
  fslmaths "s_${ICADenoised}" -bptf 31.25 -1 -add temporalmean.nii.gz -mul ${funcmask} "hpf_s_${ICADenoised}"
  done
  echo
done
