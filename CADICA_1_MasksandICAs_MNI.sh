#!/bin/sh
# CADICA - Computer-Assisted Denoising with ICA.
# It is more comprehensive than other current automatic ICA denoising in the sense it takes into account non-bold physiologic and scanner-related noise alongside motion (e.g., no need to nuissance regress WM or CSF afterwards)
# This is made to be run easily following fMRIPrep with the folder structure, but remains very feasible following any set of preprocessing, especially if generic confound parameters of interest are accessible.
# This assumes anat files and functionals are warped to mni_icbm152_nlin_asym space! This is commonly used in most pipelines now. This is different than linear or symmetric mni 152 space.
# This is modeled similarly to ALT (automated labeling tool) project, and somewhat similar to ICA-AROMA with numerous modifications.
# This uniquely takes into account different tissue types (Outbrain, GM, WM, CSF), different frequency ranges (low, BOLD, high), spikes in the time series, and significant confound correlations.
# This leads to a computer-assisted version of Manual ICA Denoising, following current gold-standard practices for IC labeling.
# CADICA can work well on its own. It also greatly increases the efficiency of Manual ICA Denoising as well as it limits the number of components an expert may want to examine for labeling.

# Decisions on Labeling ICs as Signal vs Noise is Based on the Following:
# (1) The significant clusters of an IC should significantly favor GM. Especially as compared to WM, CSF, or outside the brain.
# (2) The power of the IC should significantly favor BOLD-related physiologic frequencies (0.008-0.09 Hz). Especially as compared to lower or higher frequencies.
# (3) The IC timeseries should not show highly intense or numerous spikes (using a conservative cut off of the absolute value mean must be greater than the absolute value standard deviation)
# (4) The IC timeseries should not strongly correlate to movement confounds (all 6 motion parameters) (using a conservative cut-off of 50% of the variance of an IC accounted for by a motion confound)
# These qualifications look to help automate the current best practice protocols for the gold standard of Manual ICA Denoising.

# For the majority of data, (1) and (2) are sufficient while (3) and (4) are added safeguards that may pick up on other potentially missed noise ICs.
# There are 3 sets of selectivity with signal labeling: low, medium, and high. The user can decide based on QC what is best for their data. They can similarly decide aggressive/non-aggressive denoising.
# The labeling of signal vs noise can be easily examined and edited before applying the denoising step(s) - therefore assisting with the efficiency and accuracy of manual ICA denoising.
# The default setting is medium selectivity signal labeling (ICs labeled as signal are likely to contain true signal of interest) paired with non-aggressive denoising to retain misslabeled true signal.
# Another strong contender is low selectivity signal labeling (ICs labeled as noise are most likely to be noise) paired with aggressive denoising to remove misslabeled true noise.

# Following labeling the following is run:
# (1) non-aggressive/aggressive denoising
# (2) 6mm smoothing
# (3) Detrending (High Pass Filtering)

# Results might be improved with aggressive denoising. This is simple to change in the fsl function, or to include the CADICA noise as confounds in a GLM (e.g. in CONN Denoising)
# If doing non-aggressive denoising, consider trying to be more selective in signal labeling (high selectivity).
# If self-labeling or editing the labels (Manual Denoising), consider using the low or medium selectivity signal labels as a guide for what ICs you wish to examine yourself.
# Results have been compared to Manual Labeling, ICA-AROMA, and ALT. It appears robust, and labeling of ICs follow expectations with manual ICA denoising.
# Certain factors may perform better for your data if adjusted by the user (as discussed in notes above).
# Check for changing factors in the code to fit your data. Such as TR, filenames/locations, and scan length.

######################### Variables to Change For Your Data As Needed #############################################
# Naming of subject files is assumed to be like sub-001, sub-002, etc. (like in fMRIPrep). Therefore, you only need the numbers like 001, 002, 003 for currsubjids
# e.g. (001 002 003)
currsubjids=(108)
# currsubjids=(144)
# where your data is held. Works well with fMRIPrep folder Structure. Should be like "data/sub-###/ses-##/func/"
inputdata="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/fmriprep"
taskid="rest"
funcmaskid="rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii." # something specific in naming to the functional mask of interest
funcfileid="rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii." # something specific in naming to the functional file of interest
confoundfileid="rest_desc-confounds_timeseries." # something specific in naming to the measured confounds of the functional file of interest .tsv or .csv (ie motion parameters, global signal, csf, wm)
CADICAfuncs="/Users/keithdodd/CADICA_MNI_github" # your folder containing the CADICA scripts and such
templatefol="${CADICAfuncs}/mni_icbm152_nlin_asym_09c" # Folder containing mni templates of interest. Included in the CADICA folder
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA" # Where you want the output to be stored
cut_off_period="100" # 100 second period. This equates to 1/100 -> 0.01 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
#################################################################################################################

# Takes into account multiple sessions
for l in "${currsubjids[@]}"
do
  echo "Running for Subject ${l}"
  for j in "${sessids[@]}"
  do

  ######################### START GENERAL SET UP ##############################################################
  echo "  Running for session ${j}"
  # go to current inputdata subject folder
  currsubjfol="${inputdata}/sub-${l}"
  cd ${currsubjfol}

  # find functional folder!
  funcfol="${currsubjfol}/ses-${j}/func"

  if [ ! -d "${funcfol}" ]
  then
    echo "ERROR: Cannot find func folder holding functional files at ${funcfol}"
    echo "Exiting..."
    exit
  fi

  # go to derivatives folder & make a folder for CADICA and subject folder
  cd "${inputdata}/../"

  # make CADICA (Computer Assisted Denoising - ICA) folder if we have not already
  if [ ! -d "${CADICAfol}" ]
  then
    mkdir "${CADICAfol}"
  fi
  cd ${CADICAfol}

  # make subject folder if we have not already
  if [ ! -d "${CADICAfol}/sub-${l}" ]
  then
    cd ${CADICAfol}
    mkdir "sub-${l}"
  fi
  CADICAsubfol="${CADICAfol}/sub-${l}" # main subject folder to work within
  cd ${CADICAsubfol} # go into subject folder

  # make session folder if it does not exist
  if [ ! -d "${CADICAsubfol}/ses-${j}" ]
  then
    mkdir "${CADICAsubfol}/ses-${j}"
  fi

  sessiondir="${CADICAsubfol}/ses-${j}"
  cd ${sessiondir}

  # make session anat mask folder if we have not already
  if [ ! -d "${CADICAsubfol}/ses-${j}/anatmasks" ]
  then
    mkdir "anatmasks" # make a mask folder to put all anatomical masks of interest into
  fi
  anatmaskdir="${sessiondir}/anatmasks"

  # copy over session functional files and confounds file
  find ${funcfol}/ -name "*${funcmaskid}*" -exec cp '{}' ${sessiondir}/funcmask.nii.gz \;
  find ${funcfol}/ -name "*${funcfileid}*" -exec cp '{}' ${sessiondir}/funcfile.nii.gz \;
  find ${funcfol}/ -name "*${confoundfileid}*" -exec cp '{}' ${sessiondir}/confounds_timeseries.csv \;

  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  funcfile="${sessiondir}/funcfile.nii.gz" # functional file

  # point to anatomical aspects of interest from mni
  GMprob="${templatefol}/mni_icbm152_gm_tal_nlin_asym_09c.nii.gz"
  WMprob="${templatefol}/mni_icbm152_wm_tal_nlin_asym_09c.nii.gz"
  CSFprob="${templatefol}/mni_icbm152_csf_tal_nlin_asym_09c.nii.gz"
  anatmask="${templatefol}/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz"

  chosen_funcfile="${sessiondir}/funcfile.nii.gz"
  chosen_funcmask="${sessiondir}/funcmask.nii.gz"
  #######################################################################################################


  ##################### CREATE ROI MASKS ################################################################
  echo "  Calculating Masks in Functional Space"

  ### Edgemask: any functional outside of eroded anatomical
  # Edgemask is the only mask based on the subjects brain, the rest are based in MNI for more flexibility
  # Erode anatomical mask
  fslmaths ${anatmask} -ero "${anatmaskdir}/anatmask_eroded.nii.gz"
  # resample to functional
  flirt -ref "${chosen_funcmask}" -in "${anatmaskdir}/anatmask_eroded.nii.gz" -out "${anatmaskdir}/anatmask_eroded_resam.nii.gz" -usesqform -applyxfm
  # any functional signal outside of the eroded anat image must be edge
  fslmaths ${chosen_funcmask} -sub "${anatmaskdir}/anatmask_eroded_resam.nii.gz" -thr 0.5 "${anatmaskdir}/Edgemask_prop_final.nii.gz"

  ### GM, WM, and CSF Masks: We will calculate relative probability between the three for each voxel outside of Edgemask, based on MNI
  # Note: CSF mask really also includes spaces that would commonly have vessels as well.
  # resample each to functional space
  flirt -ref "${chosen_funcmask}" -in "${GMprob}" -out "${anatmaskdir}/GMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${chosen_funcmask}" -in "${WMprob}" -out "${anatmaskdir}/WMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${chosen_funcmask}" -in "${CSFprob}" -out "${anatmaskdir}/CSFprob_resam.nii.gz" -usesqform -applyxfm

  # Calculate relative proportion between GM, WM, and CSF, while removing Edge
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -add "${anatmaskdir}/WMprob_resam.nii.gz" -add "${anatmaskdir}/CSFprob_resam.nii.gz" "${anatmaskdir}/summedprobs_resam.nii.gz"
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -div "${anatmaskdir}/summedprobs_resam.nii.gz" -sub "${anatmaskdir}/Edgemask_prop_final.nii.gz" -thr 0.5 -mul "${chosen_funcmask}" "${anatmaskdir}/GM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -div "${anatmaskdir}/summedprobs_resam.nii.gz" -sub "${anatmaskdir}/Edgemask_prop_final.nii.gz" -thr 0.5  -mul "${chosen_funcmask}" "${anatmaskdir}/WM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/CSFprob_resam.nii.gz" -div "${anatmaskdir}/summedprobs_resam.nii.gz" -sub "${anatmaskdir}/Edgemask_prop_final.nii.gz" -thr 0.5 -mul "${chosen_funcmask}" "${anatmaskdir}/CSF_prop_final.nii.gz"

  ### WMCSF boundary, Outbrain, and Susceptibilitymask_prop_final
  # dilate and multiply CSF and WM to get the general boundary
  fslmaths "${anatmaskdir}/CSF_prop_final.nii.gz" -dilM -mul "${anatmaskdir}/WM_prop_final.nii.gz" -dilM -thr 0.5 "${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"
  # CSF plus Edgemask is all of the "outside of brain tissue"
  fslmaths "${anatmaskdir}/Edgemask_prop_final.nii.gz" -max "${anatmaskdir}/CSF_prop_final.nii.gz" "${anatmaskdir}/Outbrain_prop_final.nii.gz"
  # anatomical that does not show well in functional is susceptibility
  fslmaths "${anatmaskdir}/anatmask_resam.nii.gz" -sub "${chosen_funcmask}" -thr 0 -bin -dilF "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"

  echo "  Functional Space Masks are Computed"

  # relabel these for easier future coding
  Anatmask_resam="${anatmaskdir}/anatmask_eroded_resam.nii.gz"
  Edgemask="${anatmaskdir}/Edgemask_prop_final.nii.gz"
  GMmask="${anatmaskdir}/GM_prop_final.nii.gz"
  WMmask="${anatmaskdir}/WM_prop_final.nii.gz"
  WMCSFmask="${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"
  CSFmask="${anatmaskdir}/CSF_prop_final.nii.gz"
  Outbrainmask="${anatmaskdir}/Outbrain_prop_final.nii.gz"
  Susceptmask="${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"

  echo "    Copying of Files & Calculating Masks is Done! Now Running Melodic!"
  ##############################################################################################################


  ##################### MELODIC ICA DECOMPOSITION SECTION ######################################################
  # let's finally get to melodic IC decomposition
  cd ${sessiondir}

  # MELODIC SECTION
  # make session folder if we have not already
  if [ ! -d "${sessiondir}/melodic" ]
  then
    mkdir melodic_hpf_s # make melodic folder to run melodic in
  fi
  melfol="${sessiondir}/melodic"
  # Run Melodic on next line! Default settings. Make sure TR is correct!
  # melodic --in="${chosen_funcfile}" --outdir=${melfol} --mask="${chosen_funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}

  # grab explained variance percentage in case of use in weighting the Data
  awk '{print $1}' ${melfol}/melodic_ICstats > ${sessiondir}/IC_exp_variance.txt

  # Merge the probability maps and threshold maps which may prove useful later
  probmapnames="$(ls ${melfol}/stats/probmap_* | sort -V)"
  fslmerge -t ${melfol}/ICprobabilities.nii.gz ${probmapnames}

  thresholdnames="$(ls ${melfol}/stats/thresh_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICthresh_zstat.nii.gz ${thresholdnames}

  cd ${sessiondir}

  echo "    Melodic is Complete! Now Calculating ICA Cluster Locations & Relevant Values"
  #########################################################################################

  ####################### CALCULATE REGION OVERLAP SECTION ################################
  # Full volume, get 95% probability mask, multiple by zstat
  fslmaths ${melfol}/ICprobabilities.nii.gz -thr 0.95 -bin -mul ${melfol}/ICthresh_zstat.nii.gz -abs -mul "${chosen_funcmask}" "${sessiondir}/fullvolICA_adj.nii.gz"
  fslstats -t "${sessiondir}/fullvolICA_adj.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > fullvolume_ICmean.txt
  awk '{print $2}' tmp.txt > fullvolume_ICnumvoxels.txt

  # GM voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${GMmask}" "${sessiondir}/GMICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/GMICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > GM_ICmean.txt
  awk '{print $2}' tmp.txt > GM_ICnumvoxels.txt

  # WM voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${WMmask}" "${sessiondir}/WMICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/WMICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > WM_ICmean.txt
  awk '{print $2}' tmp.txt > WM_ICnumvoxels.txt

  # CSF voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${CSFmask}" "${sessiondir}/CSFICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/CSFICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > CSF_ICmean.txt
  awk '{print $2}' tmp.txt > CSF_ICnumvoxels.txt

  # Edge voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${Edgemask}" "${sessiondir}/EdgeICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/EdgeICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > Edge_ICmean.txt
  awk '{print $2}' tmp.txt > Edge_ICnumvoxels.txt

  # Outbrain voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${Outbrainmask}" "${sessiondir}/OutbrainICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/OutbrainICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > Outbrain_ICmean.txt
  awk '{print $2}' tmp.txt > Outbrain_ICnumvoxels.txt

  # WM CSF Boundary voxels
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${WMCSFmask}" "${sessiondir}/WMCSFICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/WMCSFICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > WMCSF_ICmean.txt
  awk '{print $2}' tmp.txt > WMCSF_ICnumvoxels.txt

  # Susceptibility voxels (Subependymal)
  fslmaths "${sessiondir}/fullvolICA_adj.nii.gz" -mul "${Susceptmask}" "${sessiondir}/SusceptICA_weighted.nii.gz"
  fslstats -t "${sessiondir}/SusceptICA_weighted.nii.gz" -M -V > ${sessiondir}/tmp.txt
  awk '{print $1}' tmp.txt > Suscept_ICmean.txt
  awk '{print $2}' tmp.txt > Suscept_ICnumvoxels.txt

  rm tmp.txt
  ###########################################################################################################


  done
done
echo
echo "1_CADICA_MasksandICAs is done running!"
echo "Next Step is to Run 2_CADICA_Labeling_MNI in Matlab!"
