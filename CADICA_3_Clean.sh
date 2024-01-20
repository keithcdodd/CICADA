#!/bin/bash

# Run this after doing auto and manaul ica selection

usage(){ echo "Usage: `basename $0` -c <cadica_dir> -p <subj_num> -s <session_ID_number> -t <task_name_include_run_labels> -i <IC_selection> -l <low_frequency_cutoff> -h <high_frequency_cutoff> -b <blur_gaussian_mm>
c: CADICA directory (necessary)
p: bids ID number for participant/subject (necessary)
s: session ID number (necessary)
t: task ID for subject (needs to include run naming) (necessary)
i: IC label selection (either manual or auto) (default is manual)
l: Low frequency cutoff in HZ (default is 0.008, there will be versions with highpass or no frequency filtering regardless)
h: High frequency cutoff in HZ (default is 0.15, there will be versions with highpass or no frequency filtering regardless)
b: Gaussian blurring in mm (default 6mm)

Runs Denoising/Cleaning of functional file, following ICA labeling from previous scripts/functions.

Example: `basename $0` -c /path/cadica_MNI152NLin6Asym_res-02/ -p 102 -s 01 -t foodpics_run-01 -f nonagg -i manual -l 0.008 -h 0.15 -b 6
" 1>&2; exit 1; }

if [ $# -lt 8 -o $# -gt 18 ]; then
  echo
  echo
	usage
  echo
  echo
fi

# Initialize things if not specified
cadica_dir="" # initialize to be edited later relative to bids_data
participant_id=""
sess_id=""
task_id=""
ICselect=""
lowfreq_cutoff=""
highfreq_cutoff=""
blur=""

while getopts "c:p:s:t:i:l:h:b:" opt; do
    case "${opt}" in
      c)
    		cadica_dir=${OPTARG}
    		;;
      p)
        participant_id=${OPTARG}
        ;;
      s)
        sess_id=${OPTARG}
        ;;
      t)
        task_id=${OPTARG}
        ;;
      i)
        ICselect=${OPTARG}
        ;;
      l)
        lowfreq_cutoff=${OPTARG}
        ;;
      h)
        highfreq_cutoff=${OPTARG}
        ;;
      b)
        blur=${OPTARG}
        ;;
      *)
        usage
        ;;
    esac
done


# OK, first check that all necessary options are given
if [ "${cadica_dir}" = "" ]; then
  echo "Missing a CADICA directory input"
  usage
  exit
fi

if [ "${participant_id}" = "" ]; then
  echo "Missing a participant id input"
  usage
  exit
fi

if [ "${sess_id}" = "" ]; then
  echo "Missing a session id input"
  usage
  exit
fi

if [ "${task_id}" = "" ]; then
  echo "Missing a task id input"
  usage
  exit
fi



# Now, check that these options all exist/work
if [ ! -d "${cadica_dir}" ]; then
  echo "CADICA Directory Not Found at ${cadica_dir}"
  echo "Exiting..."
  exit
fi

cd ${cadica_dir}

cadica_dir=$(pwd)

if [ ! -d "${cadica_dir}/sub-${participant_id}" ]; then
  echo "Participant Not Found At ${cadica_dir}/sub-${participant_id}"
  echo "Exiting..."
  exit
fi

subj_dir="${cadica_dir}/sub-${participant_id}"

if [ ! -d "${subj_dir}/ses-${sess_id}" ]; then
  echo "Session Not Found At ${subj_dir}/ses-${sess_id}"
  echo "Exiting..."
  exit
fi

sess_dir="${subj_dir}/ses-${sess_id}"

if [ ! -d "${sess_dir}/${task_id}" ]; then
  echo "Task Not Found At ${sess_dir}/ses-${task_id}"
  echo "Exiting..."
  exit
fi

# Add defaults for things that are not specified
if [ "${ICselect}" = "" ]; then
  ICselect="manual"
fi

if [ "${lowfreq_cutoff}" = "" ]; then
  lowfreq_cutoff="0.008"
fi

if [ "${highfreq_cutoff}" = "" ]; then
  work_dir="0.15"
fi

if [ "${blur}" = "" ]; then
  blur="6"
fi


# let's keep logs of everything:
date
now=`date +"%Y%m%d%H%M%S"`

log_dir="${cadica_dir}/../logs" # directory to save log file
cadica_clean_logfile="${log_dir}/cadica_clean_output_sub_${participant_id}_sess_${sess_id}_task_${task_id}_${now}.log"

(
###################### ACTUAL WORK STARTS HERE #########################################################

# give folder to ICselect
if [ "${ICselect}" = "auto" ]; then
  IC_select_dir="ic_auto_selection"
  cleaned_fol_name="cleaned_auto"
else
  IC_select_dir="ic_manual_selection"
  cleaned_fol_name="cleaned_manual"
fi

# go to current fmriprep subject folder
participant_dir="${cadica_dir}/sub-${participant_id}"
session_dir="${participant_dir}/ses-${sess_id}"
task_dir="${session_dir}/${task_id}"
cd ${task_dir}

if [ -d ${cleaned_fol_name} ]
then
  rm -rf ${cleaned_fol_name}
fi
mkdir ${cleaned_fol_name}

cd ${cleaned_fol_name}

funcfile="${task_dir}/funcfile.nii.gz"
funcmask="${task_dir}/funcmask.nii.gz"
melfol="${task_dir}/melodic"
region_mask_dir="${task_dir}/region_masks"
anatprob="${region_mask_dir}/Anatprob_resam.nii.gz"
cleaned_fol="${task_dir}/${cleaned_fol_name}"

fslmaths "${anatprob}" -thr 0.5 -bin "${region_mask_dir}/anatmask_final.nii.gz"
anatmask="${region_mask_dir}/anatmask_final.nii.gz"

# Anatomical Mask Import from before
fslmaths "${region_mask_dir}/Edge_prop.nii.gz" -sub "${region_mask_dir}/GMprob_resam.nii.gz" -thr 0.95 -bin "${cleaned_fol}/Edge_selective_mask.nii.gz"
fslmaths "${region_mask_dir}/WMandCSF_prop_final.nii.gz" -sub "${region_mask_dir}/GMprob_resam.nii.gz" -thr 0.95 -bin "${cleaned_fol}/WMCSF_selective_mask.nii.gz"
fslmaths "${funcmask}" -sub "${anatprob}" -sub "${region_mask_dir}/Edge_prop.nii.gz" -thr 0.95 -bin "${cleaned_fol}/Outbrain_noEdge_selective_mask.nii.gz"
fslmaths "${region_mask_dir}/WMprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleaned_fol}/WM_selective_mask.nii.gz"
fslmaths "${region_mask_dir}/CSFprob_resam.nii.gz" -thr 0.95 -mul "${funcmask}" -bin "${cleaned_fol}/CSF_selective_mask.nii.gz"

Edgeselectivemask="${cleaned_fol}/Edge_selective_mask.nii.gz"
WMCSFselectivemask="${cleaned_fol}/WMCSF_selective_mask.nii.gz"
Outbrainselectivemask="${cleaned_fol}/Outbrain_noEdge_selective_mask.nii.gz"
WMselectivemask="${cleaned_fol}/WM_selective_mask.nii.gz"
CSFselectivemask="${cleaned_fol}/CSF_selective_mask.nii.gz"

Edgemask="${region_mask_dir}/Edge_prop_final.nii.gz"
Subepemask="${region_mask_dir}/Subepe_prop_final.nii.gz"
CSFmask="${region_mask_dir}/CSF_prop_final.nii.gz"
OutbrainNoEdgemask="${region_mask_dir}/Outbrain_noEdge_prop_final.nii.gz"

prefix="sub-${participant_id}_ses-${sess_id}_task-${task_id}"
suffix="bold.nii.gz"

echo
nonagg_tag="CADICA_${ICselect}_nonagg"
ICADenoised_nonagg="${prefix}_${nonagg_tag}_${suffix}"
# then you run fsl_regfilt for nonaggressive or aggressive denoising! You can select your level of selection as well. This is all set at the beginning of the script.
echo "    Performing Filtering with fsl_regfilt with NonAggressive Denoising."
echo "Running: fsl_regfilt -i ${funcfile} -f $(cat ${task_dir}/${IC_select_dir}/${ICselect}_noise_dist_ICs.csv) \
-d ${melfol}/melodic_mix -m ${funcmask} -o ${ICADenoised_nonagg}"

fsl_regfilt -i ${funcfile} -f "$(cat ${task_dir}/${IC_select_dir}/${ICselect}_noise_dist_ICs.csv)" \
-d ${melfol}/melodic_mix -m ${funcmask} -o ${ICADenoised_nonagg}

#fslmaths "sub-${participant_id}_ses-${sess_id}_task-${task_id}_CADICA_${ICselect}_nonagg_bold.nii.gz" -mul ${anatmask} "sub-${participant_id}_ses-${sess_id}_task-${task_id}_CADICA_${ICselect}_nonagg_anatmasked_bold.nii.gz"


echo

echo "    Performing Filtering with Aggressive Denoising."
#echo "  Running: fsl_regfilt -i ${funcfile} -f $(cat ${task_dir}/${IC_select_dir}/${ICselect}_noise_dist_ICs.csv) \
#-d ${melfol}/melodic_mix -m ${funcmask} -a -o sub-${participant_id}_ses-${sess_id}_task-${task_id}_CADICA_${ICselect}_agg_bold.nii.gz"

#fsl_regfilt -i ${funcfile} -f "$(cat ${task_dir}/${IC_select_dir}/${ICselect}_noise_dist_ICs.csv)" \
#-d ${melfol}/melodic_mix -m ${funcmask} -a -o "sub-${participant_id}_ses-${sess_id}_task-${task_id}_CADICA_${ICselect}_agg_bold.nii.gz"

#fslmaths "sub-${participant_id}_ses-${sess_id}_task-${task_id}_CADICA_${ICselect}_agg_bold.nii.gz" -mul ${anatmask} "sub-${participant_id}_ses-${sess_id}_task-${task_id}_space-MNI152NLin2009cAsym_desc-CADICA_${ICselect}_agg_bold.nii.gz"
agg_tag="CADICA_${ICselect}_agg"
ICADenoised_agg="${prefix}_${agg_tag}_${suffix}"

# it is likely it is better, when doing aggressive denoising, to run the regression in afni (instead of fsl_regfil) with the bandpass altogether (since that is also through regression)
# Create a smoothed, highpass, and bandpass versions
3dTproject -quiet -input "${funcfile}" -ort "${task_dir}/ic_${ICselect}_selection/${ICselect}_noise_dist_covariates.txt" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -overwrite -prefix "${prefix}_${agg_tag}_${suffix}"
3dTproject -quiet -input "${funcfile}" -ort "${task_dir}/ic_${ICselect}_selection/${ICselect}_noise_dist_covariates.txt" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_${agg_tag}_${suffix}"
3dTproject -quiet -input "${funcfile}" -ort "${task_dir}/ic_${ICselect}_selection/${ICselect}_noise_dist_covariates.txt" -polort 2 -stopband 0 "${lowfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_hp_${agg_tag}_${suffix}"
3dTproject -quiet -input "${funcfile}" -ort "${task_dir}/ic_${ICselect}_selection/${ICselect}_noise_dist_covariates.txt" -polort 2 -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_bp_${agg_tag}_${suffix}"
#3dTproject -quiet -input "${ICADenoised_agg}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "s_${ICADenoised_agg}"
#3dTproject -quiet -input "${ICADenoised_agg}" -stopband 0 "${lowfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "s_hp_${ICADenoised_agg}"
#3dTproject -quiet -input "${ICADenoised_agg}" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "s_bp_${ICADenoised_agg}"

# Now do the same three steps for the nonaggressive ICA denoised files
3dTproject -quiet -input "${ICADenoised_nonagg}" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_${nonagg_tag}_${suffix}"
3dTproject -quiet -input "${ICADenoised_nonagg}" -polort 2 -stopband 0 "${lowfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_hp_${nonagg_tag}_${suffix}"
3dTproject -quiet -input "${ICADenoised_nonagg}" -polort 2 -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_bp_${nonagg_tag}_${suffix}"
echo

# calculate timeseries for CSF, WM, Edge, and just low GM in general for both aggressive and nonaggressive ICA Denoised
fslmaths "${region_mask_dir}/GMprob_resam.nii.gz" -uthr 0.05 -bin -mul ${funcmask} "${region_mask_dir}/notGM_prop_final.nii.gz"
fslmeants -i ${ICADenoised_agg} -o CSF_CADICA_agg_timeseries.txt -m ${CSFselectivemask}
fslmeants -i ${ICADenoised_agg} -o WM_CADICA_agg_timeseries.txt -m ${WMselectivemask}
fslmeants -i ${ICADenoised_agg} -o Edge_CADICA_agg_timeseries.txt -m ${Edgeselectivemask}
fslmeants -i ${ICADenoised_agg} -o notGM_CADICA_agg_timeseries.txt -m "${region_mask_dir}/notGM_prop_final.nii.gz"
fslmeants -i ${ICADenoised_agg} -o Global_CADICA_agg_timeseries.txt -m ${funcmask}

fslmeants -i ${ICADenoised_nonagg} -o CSF_CADICA_nonagg_timeseries.txt -m ${CSFselectivemask}
fslmeants -i ${ICADenoised_nonagg} -o WM_CADICA_nonagg_timeseries.txt -m ${WMselectivemask}
fslmeants -i ${ICADenoised_nonagg} -o Edge_CADICA_nonagg_timeseries.txt -m ${Edgeselectivemask}
fslmeants -i ${ICADenoised_nonagg} -o notGM_CADICA_nonagg_timeseries.txt -m "${region_mask_dir}/notGM_prop_final.nii.gz"
fslmeants -i ${ICADenoised_nonagg} -o Global_CADICA_nonagg_timeseries.txt -m ${funcmask}

fslmeants -i "${funcfile}" -o EdgeOnly_timeseries.txt -m "${Edgeselectivemask}"
fslmeants -i "${funcfile}" -o SubepeOnly_timeseries.txt -m "${WMCSFselectivemask}"
fslmeants -i "${funcfile}" -o CSFOnly_timeseries.txt -m "${CSFselectivemask}"
fslmeants -i "${funcfile}" -o OutbrainOnly_timeseries.txt -m "${Outbrainselectivemask}"

# Now combine ICA Denoised regressor files (agg and nonagg) to make design matrix in text form - in case you later want to add on more regression after ICA (like in ICA-AROMA)
paste Edge_CADICA_agg_timeseries.txt WM_CADICA_agg_timeseries.txt CSF_CADICA_agg_timeseries.txt > designregressorsICA_agg.txt
paste Edge_CADICA_nonagg_timeseries.txt WM_CADICA_nonagg_timeseries.txt CSF_CADICA_nonagg_timeseries.txt > designregressorsICA_nonagg.txt

# Also combine the funcfile regressors: You can consider just regressing these as a 4P regression of the funcfile, good in theory, but not yet well-tested
paste EdgeOnly_timeseries.txt SubepeOnly_timeseries.txt CSFOnly_timeseries.txt OutbrainOnly_timeseries.txt > designregressors.txt

echo

# Create smoothed, high passed, and bandpassed options of the original for comparisons, also copy over the original functional file for later comparison if desired
cp "${funcfile}" "./${prefix}_orig_${suffix}"
funcfile="${prefix}_orig_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_orig_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -stopband 0 "${lowfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_hp_orig_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_bp_orig_${suffix}"

# Calculate 9P regression: smoothed, highpassed, and bandpassed (pretty standard)
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/9p_regressors.txt" -mask "${funcmask}" -overwrite -prefix "${prefix}_9p_${suffix}"
3dTproject -quiet -input "${prefix}_9p_${suffix}" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_9p_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/9p_regressors.txt" -stopband 0 "${lowfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_hp_9p_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/9p_regressors.txt" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_bp_9p_${suffix}"

echo
# Calculate 8P regression: also standard, no global signal regression, smoothed, highpassed, and bandpassed
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/8p_regressors.txt" -mask "${funcmask}" -overwrite -prefix "${prefix}_8p_${suffix}"
3dTproject -quiet -input "${prefix}_8p_${suffix}" -polort 2 -mask "${task_dir}/funcmask.nii.gz" -blur "${blur}" -overwrite -prefix "${prefix}_s_8p_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/8p_regressors.txt" -stopband 0 "${lowfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_hp_8p_${suffix}"
3dTproject -quiet -input "${funcfile}" -polort 2 -ort "${task_dir}/regressors_timeseries/8p_regressors.txt" -passband "${lowfreq_cutoff}" "${highfreq_cutoff}" -mask "${funcmask}" -blur "${blur}" -overwrite -prefix "${prefix}_s_bp_8p_${suffix}"

echo

# create a masks, timeseries, regressors Folder, move all files into there, then move the bold final images back up
cd ${cleaned_fol}
mkdir masks_ts_regressors
mv *.* "${cleaned_fol}/masks_ts_regressors/"

cd "${cleaned_fol}/masks_ts_regressors/"
mv *bold.nii.g* "${cleaned_fol}/"
cd ${cleaned_fol}

echo

# Now we have fully preprocessed functionals and now can move on to analyses. For example, seed-based, or dual regression of group ICA...
) 2>&1 | tee ${cadica_clean_logfile}
