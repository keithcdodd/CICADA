#!/bin/sh
# CADICA - Complete AutoDenoise - ICA: Can follow fMRIPrep.
# "Complete" in the sense it takes into account non-bold physiologic and scanner-related noise alongside motion (no need to nuissance regress WM or CSF afterwards like ICA-AROMA)
# This is made to be run easily following fMRIPrep.
# This is modeled similar to ALT (automated labeling tool) project,
# This uniquely takes advantage of subject specific anatomy masks of WM, CSF, and Edge masks in addition to GM
# Uses FSL and Matlab.

# Decisions on Labelling ICs as Signal vs Noise is Based on the Following:
# (1) The significant clusters of an IC should lie more in Grey Matter  than WM, CSF, or outside the brain_mask
# (2) The power of the IC should lie more in BOLD-related physiologic frequencies (0.01-0.1 Hz) than below or above
# (3) The IC timeseries should not show highly intense or numerous spikes (conservative cut off)
# (4) The IC timeseries should not strongly correlate to confounds (global signal, WM, CSF, and all 6 motion parameters) (conservative cut-off)

# For the majority of data, (1) and (2) are sufficient while (3) and (4) are added safeguards that may pick up on other obvious missed noise
# For more aggressive selection, (1) and (2) can be modified from requiring "more of the signal" to requiring "a majority of the signal."
# but this might be too aggressive for particularly noisy data sets

# Following Labeling the following is run:
# (1) non-aggressive denoising
# (2) 6mm smoothing
# (3) Detrending (High Pass Filtering)

# Results might be improved with aggressive denoising. This is simple to change in the fsl function, or to include the CADICA noise as confounds in a GLM (e.g. in Conn Denoising)
# If running denoising elsewhere (e.g. CONN) we suggest also controlling for the task effect to smooth out settling in the beginning. Otherwise, consider deleting the first few volumes.
# If doing non-aggressive denoising, consider trying to be more selective in signal labeling with majority instead of simply more than (see comments within code)
# Results have been compared to Manual Labeling, ICA-AROMA, and ALT. It appears robust.
# Certain factors may perform better if adjusted by the user (as discussed in notes above).
# Check for changing factors in the code to fit your data. Such as TR, filenames/locations, and scan length.

######################### Variables to Change For Your Data As Needed #############################################
# curr subjects can be written as just the numbers as output in fMRIPrep. E.g. sub-102 is simply 102
# e.g. (001 002 003)
currsubjids=(102 103)
# where your fMRIPrep data is held
fmriprepdata="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/fmriprep"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA" # your folder containing the CADICA scripts and such
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds
#################################################################################################################

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
  # this might need slight edits depending on your naming conventions
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
      echo "Cannot find relevant anatomy folder. Compare naming structure of this code compared to your files/directories."
      exit 1
  fi

  # find functional folder too!
  funcfol="${currsubjfol}/ses-${j}/func"

  # go to derivatives folder & make a folder for CADICA and subject folder
  cd "${fmriprepdata}/../"
  derivatives="$(pwd)"

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
  cp "${funcfol}/sub-${l}_ses-${j}_task-${taskid}_desc-confounds_timeseries.tsv" ${sessiondir}/confounds_timeseries.csv

  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  funcfile="${sessiondir}/funcfile.nii.gz" # functional file

  # copy over anatomical stuff
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz" ${anatmaskdir}/anatmaskorig.nii.gz # 1 is GM, 2 is WM, 3 is CSF
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz" ${anatmaskdir}/GMprobseg.nii.gz
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz" ${anatmaskdir}/WMprobseg.nii.gz
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz" ${anatmaskdir}/CSFprobseg.nii.gz
  # also copy over the original anatomy file to main folder, just for easier reference later as desired (e.g., in CONN)
  cp "${anatfol}/${anatfilestart}_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz" ${sessiondir}/anatfile.nii.gz

  # relabel it all for easier coding
  GMprob="${anatmaskdir}/GMprobseg.nii.gz"
  WMprob="${anatmaskdir}/WMprobseg.nii.gz"
  CSFprob="${anatmaskdir}/CSFprobseg.nii.gz"
  anatmask="${anatmaskdir}/anatmaskorig.nii.gz" # mask of anat file

  # calculate in the brain and out of the brain based on the probability segments
  fslmaths "${GMprob}" -add "${WMprob}" -add "${CSFprob}" "${anatmaskdir}/inbrainprob.nii.gz"
  # for out of brain, we need to resample funcmask to same space as others for now
  flirt -ref ${anatmask} -in ${funcmask} -out "${anatmaskdir}/funcmask_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  fslmaths "${anatmaskdir}/funcmask_resam.nii.gz" -add "${anatmask}" -bin -sub "${anatmaskdir}/inbrainprob.nii.gz" -thr 0.01 "${anatmaskdir}/outbrainprob.nii.gz"

  # relabel these for coding purposes
  inbrainprob="${anatmaskdir}/inbrainprob.nii.gz"
  outbrainprob="${anatmaskdir}/outbrainprob.nii.gz"

  # resample anatseg, anatmask, anatprobs so it is all in the same space. Nearest neighbor for binary masks, trilinear for probabilities.
  flirt -ref ${funcmask} -in ${anatmask} -out "${anatmaskdir}/anatmask_final.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm -interp "nearestneighbour"
  flirt -ref ${funcmask} -in ${GMprob} -out "${anatmaskdir}/GMprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${WMprob} -out "${anatmaskdir}/WMprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${CSFprob} -out "${anatmaskdir}/CSFprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${inbrainprob} -out "${anatmaskdir}/inbrainprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm
  flirt -ref ${funcmask} -in ${outbrainprob} -out "${anatmaskdir}/outbrainprobseg_resam.nii.gz" -init "${CADICAfuncs}/eye.mat" -applyxfm

  # trilinear interp can really extend probabilities out further than is reasonable. Easiest to just threshold out the low values. 0.01 is very conservative.
  fslmaths "${anatmaskdir}/GMprobseg_resam.nii.gz" -thr 0.01 "${anatmaskdir}/GMprobseg_final.nii.gz"
  fslmaths "${anatmaskdir}/WMprobseg_resam.nii.gz" -thr 0.01 "${anatmaskdir}/WMprobseg_final.nii.gz"
  fslmaths "${anatmaskdir}/CSFprobseg_resam.nii.gz" -thr 0.01 "${anatmaskdir}/CSFprobseg_final.nii.gz"
  fslmaths "${anatmaskdir}/inbrainprobseg_resam.nii.gz" -thr 0.01 "${anatmaskdir}/inbrainprobseg_final.nii.gz"
  fslmaths "${anatmaskdir}/outbrainprobseg_resam.nii.gz" -mul "${funcmask}" -thr 0.01 "${anatmaskdir}/outbrainprobseg_final.nii.gz"

  # relabel these for easier coding
  anatbrainmask="${anatmaskdir}/anatmask_final.nii.gz" # from anat mask originally resampled
  GMprob_resam="${anatmaskdir}/GMprobseg_final.nii.gz"
  WMprob_resam="${anatmaskdir}/WMprobseg_final.nii.gz"
  CSFprob_resam="${anatmaskdir}/CSFprobseg_final.nii.gz"
  inbrainprob_resam="${anatmaskdir}/inbrainprobseg_final.nii.gz"
  outbrainprob_resam="${anatmaskdir}/outbrainprobseg_final.nii.gz"

  echo "Copying of Files & Calculating Masks is Done! Now Running Melodic!"
  # let's finally get to melodic IC decomposition
  cd ${sessiondir}
  # make session folder if we have not already
  if [ ! -d "${sessiondir}/melodic" ]
  then
    mkdir melodic # make melodic folder to run melodic in
  fi
  melfol="${sessiondir}/melodic"

  # Run Melodic! Default settings. Make sure TR is correct!
  melodic --in="${funcfile}" --outdir=${melfol} --mask=${funcmask} --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}

  cd ${sessiondir}
  # For each ICA: smooth with [3x3x3] & threshold standard (gives clusters), multiply by masks, and get useful outputs (e.g. numvoxels, mean)

  echo "Melodic is Complete! Now Calculating ICA Cluster Locations & Relevant Values"
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

  echo "1_CADICA_Masksand ICAs is done running!"
  echo "Next Step is to Run 2_CADICA_Labeling in Matlab!"

  # After this is run as expected, it is time to run the matlab file!

  # Then you run the matlab file
  # might need to add the path to Matlab ahead of time. Have not found an easy and consistent way to do this yet in bash.
  cd ${sessiondir}
  #sessiondir_ml=\'${sessiondir}\'
  /Applications/MATLAB_R2022a.app/bin/matlab -batch '2_CADICA_Labeling'
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
