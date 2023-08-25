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
currsubjids=(103)
# currsubjids=(144)
# where your data is held. Works well with fMRIPrep folder Structure. Should be like "data/sub-###/ses-##/func/"
inputdata="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/fmriprep"
taskid="rest"
funcmaskid="rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii." # something specific in naming to the functional mask of interest
funcfileid="rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii." # something specific in naming to the functional file of interest
GM_anatfileid="space-MNI152NLin2009cAsym_label-GM_probseg.nii." # something specific in naming to the GM anatfile of interest
WM_anatfileid="space-MNI152NLin2009cAsym_label-WM_probseg.nii." # something specific in naming to the WM anatfile of interest
CSF_anatfileid="space-MNI152NLin2009cAsym_label-CSF_probseg.nii." # something specific in naming to the CSF anatfile of interest
anatmaskid="space-MNI152NLin2009cAsym_desc-brain_mask.nii." # something specific in naming to the anatmask of interest
confoundfileid="rest_desc-confounds_timeseries.tsv" # something specific in naming to the measured confounds of the functional file of interest .tsv or .csv (ie motion parameters, global signal, csf, wm)
CADICAfuncs="/Users/keithdodd/CADICA_MNI_github" # your folder containing the CADICA scripts and such
templatefol="${CADICAfuncs}/mni_icbm152_nlin_asym_09c" # Folder containing mni templates of interest. Included in the CADICA folder
sessids=(01) # session numbering 01 02
TR="2" # repetition time of scan in seconds
CADICAfol="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated" # Where you want the output to be stored
low_cutoff_period="125" # 100 second period. This equates to 1/100 -> 0.008 Hz cut off frequency. Do period instead of HZ to help bash math which is limited.
lowfreq_cutoff="0.008" # in Hz
highfreq_cutoff="0.15" # in Hz
blur="6" # mm gaussian blur
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
  cd "${currsubjfol}"

  # find functional folder!
  funcfol="${currsubjfol}/ses-${j}/func"

  if [ ! -d "${funcfol}" ]
  then
    echo "ERROR: Cannot find func folder holding functional files at ${funcfol}"
    echo "Exiting..."
    exit
  fi

  # find anat folder (which could be in 1 of 2 places)
  anatfol="${currsubjfol}/anat"
  if [ ! -d "${anatfol}" ]
  then
    echo "    Anat folder not found at ${anatfol}. Checking elsewhere..."
    anatfol="${currsubjfol}/ses-${j}/anat"
  fi

  if [ ! -d "${anatfol}" ]
  then
    echo "ERROR: Cannot find anat folder holding anatomy files at ${anatfol}"
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

  # Create a bandpassed and smoothed version of the functional file, in case you want to use that
  fslmaths "${sessiondir}/funcfile.nii.gz" -Tmean temporalmean.nii.gz
  3dTproject -input "${sessiondir}/funcfile.nii.gz" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${sessiondir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "s_bp_funcfile.nii.gz"
  fslmaths "s_bp_funcfile.nii.gz" -add temporalmean.nii.gz "s_bp_funcfile.nii.gz"

  funcmask="${sessiondir}/funcmask.nii.gz" # func brain mask
  # funcfile="${sessiondir}/funcfile.nii.gz" # functional file
  funcfile="${sessiondir}/s_bp_funcfile.nii.gz" # functional file

  # copy over session anatomical CSF, WM, GM masks
  find ${anatfol}/ -name "*${GM_anatfileid}*" -exec cp '{}' ${anatmaskdir}/GM_probseg.nii.gz \;
  find ${anatfol}/ -name "*${WM_anatfileid}*" -exec cp '{}' ${anatmaskdir}/WM_probseg.nii.gz \;
  find ${anatfol}/ -name "*${CSF_anatfileid}*" -exec cp '{}' ${anatmaskdir}/CSF_probseg.nii.gz \;
  find ${anatfol}/ -name "*${anatmaskid}*" -exec cp '{}' ${anatmaskdir}/anatmask.nii.gz \;

  # point to anatomical aspects from subject space
  GMprob="${anatmaskdir}/GM_probseg.nii.gz"
  WMprob="${anatmaskdir}/WM_probseg.nii.gz"
  CSFprob="${anatmaskdir}/CSF_probseg.nii.gz"
  anatmask="${anatmaskdir}/anatmask.nii.gz"
  #######################################################################################################


  ##################### CREATE ROI MASKS ################################################################
  echo "  Calculating Masks in Functional Space"

  ### GM, WM, and CSF Masks:
  # Note: CSF mask really also includes spaces that would commonly have vessels as well.
  flirt -ref "${funcmask}" -in "${GMprob}" -out "${anatmaskdir}/GMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${funcmask}" -in "${WMprob}" -out "${anatmaskdir}/WMprob_resam.nii.gz" -usesqform -applyxfm
  flirt -ref "${funcmask}" -in "${CSFprob}" -out "${anatmaskdir}/CSFprob_resam.nii.gz" -usesqform -applyxfm

  ### Calculate Susceptibility Mask
  # anatomical that does not show well in functional (did not survive functional mask) is susceptibility, dilate if mostly outside (modal)
  flirt -ref "${funcmask}" -in "${anatmask}" -out "${anatmaskdir}/anatmask_resam.nii.gz" -usesqform -applyxfm
  # modal dilation to get susceptibilitymask, so that it is not too aggressive
  fslmaths "${anatmaskdir}/anatmask_resam.nii.gz" -sub "${funcmask}" -thr 0 -dilM -mul ${funcmask} "${anatmaskdir}/Susceptibilityprob_resam.nii.gz"
  fslmaths "${anatmaskdir}/Susceptibilityprob_resam.nii.gz" -thr 0.5 -bin "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"

  # Calculate an Inbrain mask with these probabilities
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -add "${anatmaskdir}/WMprob_resam.nii.gz" -add "${anatmaskdir}/CSFprob_resam.nii.gz" "${anatmaskdir}/Anatprob_resam.nii.gz"
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -add "${anatmaskdir}/WMprob_resam.nii.gz" "${anatmaskdir}/Inbrainprob_resam.nii.gz"


  ### Edgemask: Perimeter of the functional, need to be more selective to focus on edge
  fslmaths ${funcmask} -ero -fmean "${anatmaskdir}/eroded_smoothed_func.nii.gz"
  fslmaths ${funcmask} -sub "${anatmaskdir}/eroded_smoothed_func.nii.gz" -thr 0 "${anatmaskdir}/Edge_prop.nii.gz"
  fslmaths "${anatmaskdir}/Edge_prop.nii.gz" -thr 0.67 -bin "${anatmaskdir}/Edge_prop_final.nii.gz"
  # fslmaths "${anatmaskdir}/Edgeprob_resam.nii.gz" -max "${anatmaskdir}/funcperimeter.nii.gz" -thr 0.5 -mul "${funcmask}" -bin "${anatmaskdir}/Edge_prop_final.nii.gz"

  ### Calculate GM, WM, and CSF by subtracting out Edge (perimeter)
  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -sub "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" -bin "${anatmaskdir}/GM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -sub "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" -bin "${anatmaskdir}/WM_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/CSFprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -sub "${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz" -thr 0.5 -mul "${funcmask}" -bin "${anatmaskdir}/CSF_prop_final.nii.gz"

  ### WMCSF, GMCSF, GMWM boundaries
  # We can multiply WM and CSF and rescale it to get WMCSF boundary for Subependymal and dilate once into the WM region
  fslmaths "${anatmaskdir}/WMprob_resam.nii.gz" -mul "${anatmaskdir}/CSFprob_resam.nii.gz" "${anatmaskdir}/WMCSF_boundary_init.nii.gz"
  fslstats "${anatmaskdir}/WMCSF_boundary_init.nii.gz" -R > ${anatmaskdir}/range.txt
  minimumvalue=$(awk '{print $1}' ${anatmaskdir}/range.txt)
  maximumvalue=$(awk '{print $2}' ${anatmaskdir}/range.txt)
  fslmaths "${anatmaskdir}/WMCSF_boundary_init.nii.gz" -div "${maximumvalue}" -thr 0 "${anatmaskdir}/WMCSF_boundary_prop.nii.gz"
  fslmaths "${anatmaskdir}/WMCSF_boundary_prop.nii.gz" -thr 0.5 -bin "${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"
  rm ${anatmaskdir}/range.txt

  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -mul "${anatmaskdir}/CSFprob_resam.nii.gz" "${anatmaskdir}/GMCSF_boundary_init.nii.gz"
  fslstats "${anatmaskdir}/GMCSF_boundary_init.nii.gz" -R > ${anatmaskdir}/range.txt
  minimumvalue=$(awk '{print $1}' ${anatmaskdir}/range.txt)
  maximumvalue=$(awk '{print $2}' ${anatmaskdir}/range.txt)
  fslmaths "${anatmaskdir}/GMCSF_boundary_init.nii.gz" -div "${maximumvalue}" -thr 0.5 -bin "${anatmaskdir}/GMCSF_boundary_prop_final.nii.gz"
  rm ${anatmaskdir}/range.txt

  fslmaths "${anatmaskdir}/GMprob_resam.nii.gz" -mul "${anatmaskdir}/WMprob_resam.nii.gz" "${anatmaskdir}/GMWM_boundary_init.nii.gz"
  fslstats "${anatmaskdir}/GMWM_boundary_init.nii.gz" -R > ${anatmaskdir}/range.txt
  minimumvalue=$(awk '{print $1}' ${anatmaskdir}/range.txt)
  maximumvalue=$(awk '{print $2}' ${anatmaskdir}/range.txt)
  fslmaths "${anatmaskdir}/GMWM_boundary_init.nii.gz" -div "${maximumvalue}" -thr 0.5 -bin "${anatmaskdir}/GMWM_boundary_prop_final.nii.gz"
  rm ${anatmaskdir}/range.txt

  # Subependymal space can be the WMCSF boundary, but add on a layer of WM too
  fslmaths "${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz" -dilM -mul "${anatmaskdir}/WM_prop_final.nii.gz" -thr 0.5 -bin "${anatmaskdir}/WMCSF_WMdilation.nii.gz"
  fslmaths "${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz" -max "${anatmaskdir}/WMCSF_WMdilation.nii.gz" "${anatmaskdir}/Subepe_prop_final.nii.gz"


  ### Calculate Inbrain & Outbrain: funcmask outside of anatomical inbrain+CSF+WMCSF.
  fslmaths "${funcmask}" -sub "${anatmaskdir}/anatmask_resam.nii.gz" -thr 0 -bin -dilD -bin -sub "${anatmaskdir}/Inbrainprob_resam.nii.gz" -thr 0.67 -bin -mul "${funcmask}" "${anatmaskdir}/OutbrainOnly_prop_final.nii.gz"
  fslmaths "${anatmaskdir}/Inbrainprob_resam.nii.gz" -sub "${anatmaskdir}/Edge_prop_final.nii.gz" -thr 0.5 -bin -mul "${funcmask}" "${anatmaskdir}/Inbrain_prop_final.nii.gz"
  fslmaths "${funcmask}" -sub "${anatmaskdir}/Inbrainprob_resam.nii.gz" -thr 0.5 -bin -mul "${funcmask}" "${anatmaskdir}/Outbrain_prop_final.nii.gz"

  echo "  Functional Space Masks are Computed"
  # relabel these for easier future coding
  Anatmask_resam="${anatmaskdir}/anatmask_eroded_resam.nii.gz"
  Edgemask="${anatmaskdir}/Edge_prop_final.nii.gz"
  GMmask="${anatmaskdir}/GM_prop_final.nii.gz"
  WMmask="${anatmaskdir}/WM_prop_final.nii.gz"
  CSFmask="${anatmaskdir}/CSF_prop_final.nii.gz"
  Outbrainmask="${anatmaskdir}/Outbrain_prop_final.nii.gz"
  Susceptmask="${anatmaskdir}/Susceptibilitymask_prop_final.nii.gz"
  Inbrainmask="${anatmaskdir}/Inbrain_prop_final.nii.gz"
  Subepemask="${anatmaskdir}/Subepe_prop_final.nii.gz"
  WMCSFmask="${anatmaskdir}/WMCSF_boundary_prop_final.nii.gz"
  GMCSFmask="${anatmaskdir}/GMCSF_boundary_prop_final.nii.gz"
  GMWMmask="${anatmaskdir}/GMWM_boundary_prop_final.nii.gz"


  echo "    Copying of Files & Calculating Masks is Done! Now Running Melodic!"
  ##############################################################################################################


  ##################### MELODIC ICA DECOMPOSITION SECTION ######################################################
  # let's finally get to melodic IC decomposition
  cd ${sessiondir}

  # MELODIC SECTION
  # make session folder if we have not already
  if [ -d "${sessiondir}/melodic" ]
  then
    rm -rf melodic
  fi
  mkdir melodic # make melodic folder to run melodic in
  melfol="${sessiondir}/melodic"
  # Run Melodic on next line! Default settings. Make sure TR is correct!
  melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}

  # Merge the probability maps and threshold maps which may prove useful later
  probmapnames="$(ls ${melfol}/stats/probmap_* | sort -V)"
  fslmerge -t ${melfol}/ICprobabilities.nii.gz ${probmapnames}

  thresholdnames="$(ls ${melfol}/stats/thresh_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICthresh_zstat.nii.gz ${thresholdnames}

  cd ${sessiondir}

  echo "    Melodic is Complete! Now Calculating ICA Cluster Locations & Relevant Values"
  #########################################################################################

  ## calculate and grab the smoothness values for each zstat

  # I would create a new folder in session folder first to do all this
  clusterfol="${sessiondir}/clustering"
  if [ -d "${clusterfol}" ]
  then
    rm -rf ${clusterfol}
  fi
  mkdir "${clusterfol}"
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
  done

  # save dlh array to a textfile
  printf "%s\n" "${dlh_numbers[@]}" > ${clusterfol}/smoothness.txt
  # save clustersize array to a textfile
  printf "%s\n" "${clustersizes[@]}" > ${clusterfol}/clustersizes.txt

  # combine clusterthresh_zstats into one image, this is what we will use for mapping
  clusterthresholdnames="$(ls ${clusterfol}/clusterthresh_zstat* | sort -V)"
  fslmerge -t ${melfol}/ICclusterthresh_zstat.nii.gz ${clusterthresholdnames}

  cd ${sessiondir}
  if [ -d "ROIcalcs" ]
  then
    rm -rf ROIcalcs
  fi
  mkdir "ROIcalcs"
  ROIcalcfol="${sessiondir}/ROIcalcs"

  # grab explained variance percentage in case of use in weighting the Data
  awk '{print $1}' ${melfol}/melodic_ICstats > ${ROIcalcfol}/IC_exp_variance.txt

  # with the cluster thresholded maps
  fslmaths ${melfol}/ICprobabilities.nii.gz -thr 0.95 -bin "${ROIcalcfol}/highprobmask.nii.gz"
  fslmaths ${melfol}/ICclusterthresh_zstat.nii.gz -abs "${ROIcalcfol}/fullvolICA_adj.nii.gz"
  # calculate without clusterizing
  fslmaths ${melfol}/ICthresh_zstat.nii.gz -abs -thr 3.09 "${ROIcalcfol}/fullvolICA_noclustering.nii.gz"

  # Compare before and after clustering
  fslstats -t "${ROIcalcfol}/fullvolICA_noclustering.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_noclustering_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_noclustering_ICnumvoxels.txt

  # Calculate voxel stats for after cluster correction (so you can compare, to give you an idea of clustering of IC later)
  fslstats -t "${ROIcalcfol}/fullvolICA_adj.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/fullvolume_ICnumvoxels.txt

  calcfile="${ROIcalcfol}/fullvolICA_adj.nii.gz"

  # GM voxels
  fslmaths "${calcfile}" -mul "${GMmask}" "${ROIcalcfol}/GMICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/GMICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GM_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GM_ICnumvoxels.txt

  # WM voxels
  fslmaths "${calcfile}" -mul "${WMmask}" "${ROIcalcfol}/WMICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/WMICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WM_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WM_ICnumvoxels.txt

  # CSF voxels
  fslmaths "${calcfile}" -mul "${CSFmask}" "${ROIcalcfol}/CSFICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/CSFICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/CSF_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/CSF_ICnumvoxels.txt

  # Edge voxels
  fslmaths "${calcfile}" -mul "${Edgemask}" "${ROIcalcfol}/EdgeICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/EdgeICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Edge_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Edge_ICnumvoxels.txt

  # Outbrain voxels
  fslmaths "${calcfile}" -mul "${Outbrainmask}" "${ROIcalcfol}/OutbrainICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/OutbrainICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Outbrain_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Outbrain_ICnumvoxels.txt

  # Inbrain voxels
  fslmaths "${calcfile}" -mul "${Inbrainmask}" "${ROIcalcfol}/InbrainICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/InbrainICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Inbrain_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Inbrain_ICnumvoxels.txt

  # WM CSF Boundary + surrounding WM voxels (Subependymal)
  fslmaths "${calcfile}" -mul "${Subepemask}" "${ROIcalcfol}/SubepeICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/SubepeICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Subepe_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Subepe_ICnumvoxels.txt

  # Susceptibility voxels
  fslmaths "${calcfile}" -mul "${Susceptmask}" "${ROIcalcfol}/SusceptICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/SusceptICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Suscept_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/Suscept_ICnumvoxels.txt

  # WMCSF Boundary voxels
  fslmaths "${calcfile}" -mul "${WMCSFmask}" "${ROIcalcfol}/WMCSFICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/WMCSFICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WMCSF_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/WMCSF_ICnumvoxels.txt

  # GMCSF Boundary voxels
  fslmaths "${calcfile}" -mul "${GMCSFmask}" "${ROIcalcfol}/GMCSFICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/GMCSFICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GMCSF_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GMCSF_ICnumvoxels.txt

  # GMWM Boundary voxels
  fslmaths "${calcfile}" -mul "${GMWMmask}" "${ROIcalcfol}/GMWMICA_weighted.nii.gz"
  fslstats -t "${ROIcalcfol}/GMWMICA_weighted.nii.gz" -M -V > ${ROIcalcfol}/tmp.txt
  awk '{print $1}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GMWM_ICmean.txt
  awk '{print $2}' ${ROIcalcfol}/tmp.txt > ${ROIcalcfol}/GMWM_ICnumvoxels.txt


  rm ${ROIcalcfol}/tmp.txt
  ###########################################################################################################

  done
done
echo
echo "1_CADICA_MasksandICAs is done running!"
echo "Next Step is to Run 2_CADICA_Labeling_MNI in Matlab!"
