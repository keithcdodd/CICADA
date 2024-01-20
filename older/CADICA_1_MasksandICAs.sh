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
# (1) The significant clusters of an IC should significantly favor GM. Especially as compared to Subependymal, CSF, or outside the brain.
# (2) The power of the IC should significantly favor BOLD-related physiologic frequencies (0.008-0.09 Hz). Especially as compared to lower or higher frequencies.
# (3) The IC timeseries should not show highly intense or numerous spikes (using a conservative cut off of the absolute value mean must be greater than the absolute value standard deviation)
# (4) The IC timeseries should not strongly correlate to movement confounds (all 6 motion parameters) (using a conservative cut-off of 50% of the variance of an IC accounted for by a motion confound)
# These qualifications look to help automate the current best practice protocols for the gold standard of Manual ICA Denoising.

# Code generally follows Griffanti et al guidelines for Manual ICA Denoising. IC's that could be considered "close" to being labeled signal vs noise, is recorded in matlab script to greatly diminish
# the number of ICs one may examine to check validity for the ICA Denoising.

# Following labeling the following is run:
# (1) non-aggressive/aggressive denoising
# (2) 6mm smoothing
# (3) Detrending (High Pass Filtering)
# Note: Low pass filtering is not performed as some true BOLD GM signal may be contained within the higher frequencies.

# Results might be improved with aggressive denoising. This is simple to change in the fsl function, or to include the CADICA noise as confounds in a GLM (e.g. in CONN Denoising)
# Results have been compared to Manual Labeling, ICA-AROMA, and ALT. It appears robust, and labeling of ICs follow expectations with manual ICA denoising.
# Certain factors may perform better for your data if adjusted by the user. The largest potential changes to be made are likely in the 2nd script (matlab file)
# Check for changing factors in the code to fit your data. Such as TR, filenames/locations, and scan length.

######################### Variables to Change For Your Data As Needed #############################################
# Naming of subject files is assumed to be like sub-001, sub-002, etc. (like in fMRIPrep). Therefore, you only need the numbers like 001, 002, 003 for currsubjids
# e.g. (001 002 003)
# currsubjids=(102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187)
currsubjids=(102)
# where your data is held. Works well with fMRIPrep folder Structure. Should be like "data/sub-###/ses-##/func/"
inputdata="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/fmriprep"
taskid="rest"
funcmaskid="rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii." # something specific in naming to the functional mask of interest
funcfileid="rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii." # something specific in naming to the functional file of interest
GM_anatfileid="space-MNI152NLin2009cAsym_label-GM_probseg.nii." # something specific in naming to the GM anatfile of interest
WM_anatfileid="space-MNI152NLin2009cAsym_label-WM_probseg.nii." # something specific in naming to the WM anatfile of interest
CSF_anatfileid="space-MNI152NLin2009cAsym_label-CSF_probseg.nii." # something specific in naming to the CSF anatfile of interest
anatmaskid="space-MNI152NLin2009cAsym_desc-brain_mask.nii." # something specific in naming to the anatmask of interest
anatfileid="space-MNI152NLin2009cAsym_desc-preproc_T1w.nii." # something specific in naming to the anatfile of interest
confoundfileid="rest_desc-confounds_timeseries." # something specific in naming to the measured confounds of the functional file of interest .tsv or .csv (ie motion parameters, global signal, csf, wm)
CADICAfuncs="/Users/keithdodd/CADICA_MNI_github" # your folder containing the CADICA scripts and such
templatefol="${CADICAfuncs}/mni_icbm152_nlin_asym_09c" # Folder containing mni templates of interest. Included in the CADICA folder
sessids=(01) # session numbering 01 02
TR="2" # repetition time of scan in seconds
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated" # Where you want the output to be stored
low_cutoff_period="125" # 100 second period. This equates to 1/100 -> 0.008 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
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
    echo "    ERROR: Cannot find func folder holding functional files at ${funcfol}"
    echo "    Exiting..."
    exit
  fi

  # find anat folder (which could be in 1 of 3 places)
  anatfol="${currsubjfol}/anat"
  if [ ! -d "${anatfol}" ]
  then
    echo "     Anat folder not found at ${anatfol}. Checking elsewhere..."
    anatfol="${currsubjfol}/ses-${j}/anat"
  fi

  ###(((
  # You might not keep this if statement, this is just because my second session folder does not have an anatomical
  if [ ! -d "${anatfol}" ]
  then
    echo "    Cannot find anat folder holding anatomy files at ${anatfol}. Checking first session..."
    anatfol="${currsubjfol}/ses-01/anat"
  fi
  ###)))

  if [ ! -d "${anatfol}" ]
  then
    echo "    ERROR: Cannot find anat folder holding anatomy files at ${anatfol}."
    echo "    Exiting..."
    exit
  fi

  echo "    Anat and Func folders have been located"

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
  find ${funcfol}/ -name "*${confoundfileid}*sv" -exec cp '{}' ${sessiondir}/confounds_timeseries.csv \;

  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  funcfile="${sessiondir}/funcfile.nii.gz" # functional file

  # copy over session anatomical CSF, WM, GM masks, anatmask, and anatfile
  find ${anatfol}/ -name "*${GM_anatfileid}*" -exec cp '{}' ${anatmaskdir}/GM_probseg.nii.gz \;
  find ${anatfol}/ -name "*${WM_anatfileid}*" -exec cp '{}' ${anatmaskdir}/WM_probseg.nii.gz \;
  find ${anatfol}/ -name "*${CSF_anatfileid}*" -exec cp '{}' ${anatmaskdir}/CSF_probseg.nii.gz \;
  find ${anatfol}/ -name "*${anatmaskid}*" -exec cp '{}' ${anatmaskdir}/anatmask.nii.gz \;
  find ${anatfol}/ -name "*${anatfileid}*" -exec cp '{}' ${sessiondir}/anatfile.nii.gz \;

  # point to anatomical aspects from subject space
  GMprob="${anatmaskdir}/GM_probseg.nii.gz"
  WMprob="${anatmaskdir}/WM_probseg.nii.gz"
  CSFprob="${anatmaskdir}/CSF_probseg.nii.gz"
  anatmask="${anatmaskdir}/anatmask.nii.gz"
  #######################################################################################################


  ##################### CREATE ROI MASKS ################################################################
  echo "    Calculating Masks in Functional Space"

  ### GM, WM, and CSF Masks:
  # Note: CSF mask really also includes spaces that would commonly have vessels as well.
  flirt -ref "${funcmask}" -in "${GMprob}" -out "${anatmaskdir}/GMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${funcmask}" -in "${WMprob}" -out "${anatmaskdir}/WMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${funcmask}" -in "${CSFprob}" -out "${anatmaskdir}/CSFprob_resam.nii.gz" -usesqform -applyxfm

  ### Calculate Susceptibility Mask
  # anatomical that does not show well in functional (did not survive functional mask) is susceptibility, dilate if mostly outside (modal)
  flirt -ref "${funcmask}" -in "${anatmask}" -out "${anatmaskdir}/anatmask_resam.nii.gz" -usesqform -applyxfm
  fslmaths "${funcmask}" -dilF "${anatmaskdir}/dilated_funcmask.nii.gz"
  fslmaths "${anatmaskdir}/anatmask_resam.nii.gz" -mul ${funcmask} -thr 0 "${anatmaskdir}/anatmask_resam_final.nii.gz"
  # modal dilation to get susceptibilitymask, so that it is not too aggressive
  fslmaths "${anatmaskdir}/anatmask_resam.nii.gz" -sub "${funcmask}" -thr 0 -dilD "${anatmaskdir}/Susceptibilityprob_resam.nii.gz"
  # Need to focus on frontal sinus region though
  rm "${anatmaskdir}/frontbrainfunc.nii.gz"
  3dcalc -a "${anatmaskdir}/dilated_funcmask.nii.gz" -expr 'a*isnegative(y)' -prefix "${anatmaskdir}/frontbrainfunc.nii.gz"
  fslmaths  "${anatmaskdir}/Susceptibilityprob_resam.nii.gz" -mul "${anatmaskdir}/frontbrainfunc.nii.gz" "${anatmaskdir}/Susceptibilityprobfront_resam.nii.gz"
  fslstats "${anatmaskdir}/Susceptibilityprobfront_resam.nii.gz" -R > ${anatmaskdir}/rangesusc.txt
  minimumvalue=$(awk '{print $1}' ${anatmaskdir}/rangesusc.txt)
  maximumvalue=$(awk '{print $2}' ${anatmaskdir}/rangesusc.txt)
  fslmaths "${anatmaskdir}/Susceptibilityprobfront_resam.nii.gz" -div ${maximumvalue} -thr 0.5 -mul "${funcmask}" "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"

  # Calculate an Inbrain mask with these probabilities
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -add "${anatmaskdir}/WMprob_resam.nii.gz" -add "${anatmaskdir}/CSFprob_resam.nii.gz" "${anatmaskdir}/Anatprob_resam.nii.gz"

  ### Edgemask: any functional outside of anatomical inbrain and/or the perimeter of the functional image
  fslmaths ${funcmask} -sub "${anatmaskdir}/Anatprob_resam.nii.gz" -thr 0.5 "${anatmaskdir}/Edgeprob_resam.nii.gz"
  fslmaths ${funcmask} -ero -fmean "${anatmaskdir}/eroded_funcmask.nii.gz"
    # make threshold higher so that edge does not bleed too far into the brain since erode erodes deeply in fsl
  fslmaths ${funcmask} -sub "${anatmaskdir}/eroded_funcmask.nii.gz" -thr 0.67 "${anatmaskdir}/funcbrainperimeter.nii.gz"
  fslmaths "${anatmaskdir}/Edgeprob_resam.nii.gz" -max "${anatmaskdir}/funcbrainperimeter.nii.gz" -thr 0.5 -mul "${funcmask}" "${anatmaskdir}/Edge_prop_final.nii.gz"


  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" "${anatmaskdir}/GM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" "${anatmaskdir}/WM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/CSFprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" "${anatmaskdir}/CSF_prop_final.nii.gz"

  ### Internal veins hide within WMCSF boundary, and can extend a bit more into the WM from there
  # We can add WM and CSF probability, subtract out WM and CSF thresholds, dilate to extend WM (and remove all other region mask overlap), and rescale it to get WMCSF boundary for Subependymal
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -mul "${anatmaskdir}/CSFprob_resam.nii.gz" -thr 0 "${anatmaskdir}/WMCSF_boundary.nii.gz"
  fslstats "${anatmaskdir}/WMCSF_boundary.nii.gz" -R > ${anatmaskdir}/rangewmcsf.txt
  minimumvalue=$(awk '{print $1}' ${anatmaskdir}/rangewmcsf.txt)
  maximumvalue=$(awk '{print $2}' ${anatmaskdir}/rangewmcsf.txt)
  fslmaths "${anatmaskdir}/WMCSF_boundary.nii.gz" -div "${maximumvalue}" -mul "${funcmask}" -thr 0 "${anatmaskdir}/WMCSF_boundary_prop_interim.nii.gz"
  fslmaths "${anatmaskdir}/WMCSF_boundary_prop_interim.nii.gz" -dilF -mul "${anatmaskdir}/WM_prop_final.nii.gz" -max "${anatmaskdir}/WMCSF_boundary_prop_interim.nii.gz" -thr 0.5 "${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"


  # CSF plus Edgemask (including susceptibility) is all of the "outside of brain tissue"
  fslmaths "${anatmaskdir}/Edge_prop_final.nii.gz" -max "${anatmaskdir}/CSF_prop_final.nii.gz" "${anatmaskdir}/Outbrain_prop_final.nii.gz"
  # Inbrain is GM plus WM
  fslmaths "${anatmaskdir}/GM_prop_final.nii.gz" -max "${anatmaskdir}/WM_prop_final.nii.gz" "${anatmaskdir}/Inbrain_prop_final.nii.gz"

  echo "    Functional Space Masks are Computed"
  # relabel these for easier future coding
  Anatmask_resam="${anatmaskdir}/anatmask_eroded_resam.nii.gz"
  Edgemask="${anatmaskdir}/Edge_prop_final.nii.gz"
  GMmask="${anatmaskdir}/GM_prop_final.nii.gz"
  WMmask="${anatmaskdir}/WM_prop_final.nii.gz"
  CSFmask="${anatmaskdir}/CSF_prop_final.nii.gz"
  Outbrainmask="${anatmaskdir}/Outbrain_prop_final.nii.gz"
  Susceptmask="${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"
  Inbrainmask="${anatmaskdir}/Inbrain_prop_final.nii.gz"
  WMCSFmask="${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"


  echo "    Copying of Files & Calculating Masks is Done! Now Running Melodic!"
  ##############################################################################################################

  # Try bandpass and smoothing file before doing ICA MELODIC
  # smooth 6mm and bandpass 0.008 to 0.15 HZ (conservative options)
  # prep masks - smooth 6 mm resampled and take top 98% (fsl robust range)- make design regressors from those
  fslmaths "${anatmaskdir}/Edgeprob_resam.nii.gz" -thrP 100 -mul "${funcmask}" -bin "${anatmaskdir}/HighEdgeOnly_mask.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -thrP 100 -mul "${funcmask}" -bin "${anatmaskdir}/HighWMOnly_mask.nii.gz"
  fslmaths "${anatmaskdir}/CSFprob_resam.nii.gz" -thrP 100 -mul "${funcmask}" -bin "${anatmaskdir}/HighCSFOnly_mask.nii.gz"

  fslmeants -i ${funcfile} -o Edge_timeseries.txt -m "${anatmaskdir}/HighEdgeOnly_mask.nii.gz"
  fslmeants -i ${funcfile} -o WM_timeseries.txt -m "${anatmaskdir}/HighWMOnly_mask.nii.gz"
  fslmeants -i ${funcfile} -o CSF_timeseries.txt -m "${anatmaskdir}/HighCSFOnly_mask.nii.gz"

  paste Edge_timeseries.txt WM_timeseries.txt CSF_timeseries.txt > designregressors.1D

  # would potentially include -ort designregressors.1D in there normally

  # For all afni stuff, the file must not exist first (won't overwrite)
  rm "${sessiondir}/temporalmean.nii.gz"
  rm "${sessiondir}/bp_demean_funcfile.nii.gz"
  rm "${sessiondir}/bp_funcfile.nii.gz"
  3dTstat -mean -mask "${funcmask}" -prefix "${sessiondir}/temporalmean.nii.gz" "${funcfile}"
  3dTproject -mask "${funcmask}" -input "${funcfile}" -passband 0.008 0.15 -quiet -prefix "${sessiondir}/bp_demean_funcfile.nii.gz"
  fslmaths "${sessiondir}/bp_demean_funcfile.nii.gz" -add temporalmean.nii.gz "${sessiondir}/bp_funcfile.nii.gz"

  funcfile="${sessiondir}/bp_funcfile.nii.gz"


  ##################### MELODIC ICA DECOMPOSITION SECTION ######################################################
  # let's finally get to melodic IC decomposition
  cd ${sessiondir}

  # MELODIC SECTION
  # make session folder if we have not already
  # melodic will run if session folder does not exist (will NOT rerun)

  if [ ! -d "${sessiondir}/melodic" ]
  then
    mkdir melodic # make melodic folder to run melodic in
    melfol="${sessiondir}/melodic"
    echo "    Melodic is Running!"
    # Run Melodic on next line! Default settings. Make sure TR is correct!
    melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
  fi
  melfol="${sessiondir}/melodic"


  # Merge the probability maps and threshold maps which may prove useful later
  probmapnames="$(ls ${melfol}/stats/probmap_* | sort -V)"
  fslmerge -t ${melfol}/ICprobabilities.nii.gz ${probmapnames}

  thresholdnames="$(ls ${melfol}/stats/thresh_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICthresh_zstat.nii.gz ${thresholdnames}

  cd ${sessiondir}

  echo "    Calculating ICA Cluster Locations & Relevant Values"
  #########################################################################################

  ## calculate and grab the smoothness values for each zstat

  # I would create a new folder in session folder first to do all this
  clusterfol="${sessiondir}/clustering"
  rm -rf "${clusterfol}"
  if [ ! -d "${clusterfol}" ]
  then
    mkdir "${clusterfol}"
  fi
  cd ${clusterfol}


  # first need to know number of voxels in funcmask
  fslstats -t ${funcmask} -V > ${clusterfol}/funcmask_numvoxels.txt
  funcmask_numvoxels=$(awk '{print $1}' ${clusterfol}/funcmask_numvoxels.txt)

  # create an array from files in numerical order
  thresholdarray=($(find ${melfol}/stats -name thresh_zstat*.nii.gz | sort -V))

  # declare an empty array to store the dlh numbers and clustersizes
  dlh_numbers=()
  clustersizes=()

  # Iterate over threshold array
  for (( j=0; j<${#thresholdarray[@]}; j++ ));
  do
    # grab curr file
    file="${thresholdarray[$j]}"

    # run smoothest command and extract DLH value
    smoothest -z ${file} -m ${funcmask} > ${clusterfol}/curr_smoothest.txt

    # Extract the number following the DLH label
    dlh_number=$(awk '/^DLH/{print $NF}' ${clusterfol}/curr_smoothest.txt)

    # add the DLH number to the array
    dlh_numbers+=("$dlh_number")

    # run cluster with dlh_number and funcmask volume
    cluster -i ${file} -t 3.09 -d ${dlh_number} --volume=${funcmask_numvoxels} -p 0.05 --minclustersize --osize="${clusterfol}/sizeimage_${j}.nii.gz" --no_table > ${clusterfol}/curr_cluster.txt

    # save min cluster size
    minclustersize=$(awk '/^Minimum/{print $NF}' ${clusterfol}/curr_cluster.txt)

    # add minimum cluster size to array
    clustersizes+=("$minclustersize")

    # threshold size image by minimum cluster size
    fslmaths ${clusterfol}/sizeimage_${j}.nii.gz -thr ${minclustersize} -bin -mul ${file} "${clusterfol}/clusterthresh_zstat${j}.nii.gz"

    # Also, while here, threshold to peaks in zstat image:
    fslmaths ${file} -abs "${clusterfol}/abs_zstat${j}.nii.gz"
    mean=$(fslstats "${clusterfol}/abs_zstat${j}.nii.gz" -M)
    standdev=$(fslstats "${clusterfol}/abs_zstat${j}.nii.gz" -S)
    peakcutoff=$(echo "${mean}+2*${standdev}" | bc)
    # could do -thr ${peakcutoff}
    fslmaths "${clusterfol}/abs_zstat${j}.nii.gz" -thr ${peakcutoff} "${clusterfol}/peaked_zstat${j}.nii.gz"
  done

  # save dlh array to a textfile
  printf "%s\n" "${dlh_numbers[@]}" > ${clusterfol}/smoothness.txt
  # save clustersize array to a textfile
  printf "%s\n" "${clustersizes[@]}" > ${clusterfol}/clustersizes.txt

  # combine peaked_zstats into one image
  peakedthresholdnames="$(ls ${clusterfol}/peaked_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICpeaked_zstat.nii.gz ${peakedthresholdnames}

  # combine clusterthresholded_zstats too, this is what we will use for mapping
  clusterthresholdnames="$(ls ${clusterfol}/clusterthresh_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICclusterthresh_zstat.nii.gz ${clusterthresholdnames}

  cd ${sessiondir}
  if [ ! -d "ROIcalcs" ]
  then
    mkdir "ROIcalcs"
  fi
  ROIcalcfol="${sessiondir}/ROIcalcs"

  # grab explained variance percentage in case of use in weighting the Data
  awk '{print $1}' ${melfol}/melodic_ICstats > ${ROIcalcfol}/IC_exp_variance.txt

  # with the cluster thresholded maps
  fslmaths ${melfol}/ICprobabilities.nii.gz -thr 0.999 -bin "${ROIcalcfol}/highprobmask.nii.gz"
  fslmaths ${melfol}/ICclusterthresh_zstat.nii.gz -abs "${ROIcalcfol}/fullvolICA_adj.nii.gz"
  # calculate without clusterizing
  fslmaths ${melfol}/ICthresh_zstat.nii.gz -abs -thr 3.09 "${ROIcalcfol}/fullvolICA_noclustering.nii.gz"
  fslmaths "${ROIcalcfol}/fullvolICA_adj.nii.gz" -mul "${ROIcalcfol}/highprobmask.nii.gz" "${ROIcalcfol}/highprobfullvolICA.nii.gz"
  fslmaths ${melfol}/ICpeaked_zstat.nii.gz -bin -mul ${melfol}/ICthresh_zstat.nii.gz -abs "${ROIcalcfol}/fullvolICA_peaked.nii.gz"

  # Create a "peak locator" map. This can help distinguish between difficult to tell ICs


  # Compare before and after clustering
  fslstats -t "${ROIcalcfol}/fullvolICA_noclustering.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_noclustering_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_noclustering_ICnumvoxels.txt

  # Calculate voxel stats for after cluster correction (so you can compare, to give you an idea of clustering of IC later)
  fslstats -t "${ROIcalcfol}/fullvolICA_adj.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_ICnumvoxels.txt

  # Calculate a mask for all the masks, so calculate a relative full volume
  # fslmaths ${GMmask} -add ${Edgemask} -add ${WMmask} -add ${CSFmask} -bin

  comparisonfunc="${ROIcalcfol}/fullvolICA_adj.nii.gz"

  # GM voxels
  fslmaths "${GMmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/GMICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/GMICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GM_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GM_ICnumvoxels.txt

  # WM voxels
  fslmaths "${WMmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/WMICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/WMICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WM_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WM_ICnumvoxels.txt

  # CSF voxels
  fslmaths "${CSFmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/CSFICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/CSFICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/CSF_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/CSF_ICnumvoxels.txt

  # Edge voxels
  fslmaths "${Edgemask}" -mul "${comparisonfunc}" "${ROIcalcfol}/EdgeICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/EdgeICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Edge_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Edge_ICnumvoxels.txt

  # Outbrain voxels
  fslmaths "${Outbrainmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/OutbrainICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/OutbrainICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Outbrain_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Outbrain_ICnumvoxels.txt

  # Inbrain voxels
  fslmaths "${Inbrainmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/InbrainICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/InbrainICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Inbrain_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Inbrain_ICnumvoxels.txt

  # WM CSF Boundary voxels (Subependymal)
  fslmaths "${WMCSFmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/WMCSFICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/WMCSFICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WMCSF_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WMCSF_ICnumvoxels.txt

  # Susceptibility voxels
  fslmaths "${Susceptmask}" -mul "${comparisonfunc}" "${ROIcalcfol}/SusceptICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/SusceptICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Suscept_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Suscept_ICnumvoxels.txt

  rm ${ROIcalcfol}/tmp.txt
  ###########################################################################################################

  done
done
echo
echo "1_CADICA_MasksandICAs is done running!"
echo "Next Step is to Run 2_CADICA_Labeling_MNI in Matlab!"
