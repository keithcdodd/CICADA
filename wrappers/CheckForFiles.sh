#!/bin/bash

# Run first script/function of CADICA
# NOTE: make sure your version of FSl is creating .nii.gz files before opening Matlab and running any scripts.
# e.g.,: in bash_profile:
# FSLOUTPUTTYPE=NIFTI_GZ
# export FSLOUTPUTTYPE

# immediately exit upon common error
set -euo pipefail

# initialize things
# Initialize things if not specified
input_dir=""
participant_id=""
sess_id=""
task_id=""
tr=""
brain_network_template=""
template_tag="MNI152NLin6Asym"
output_dir=""
redoMELODIC=""
anat_sess="01"
script_dir="$(dirname $0)"
brain_network_template=""
num_ICs=""

usage(){
>&2 cat << EOF
Usage: `basename $0` --input_dir='fmriprep_dir' --participant_id='part_number' --sess_id='sess_number' --task_id='task_name_with_run' --anat_sess='best_anat_sess_number' --output_dir='cicada_dir' --brain_network_template='brain_network_template.nii.gz' --num_ICs='number_ICs_for_MELODIC' --redoMELODIC
  [ -f | --input_dir input : (necessary) directory with fmriprep data (our in the same format/naming) ]
  [ -p | --participant_id input : (necessary) bids ID for participant/subject ]
  [ -s | --sess_id input : (necessary) session ID number ]
  [ -t | --task_id input : (necessary) task ID for subject (needs to include run naming) ]
  [ -l | --template_tag input : (necessary) Template Name that you are using (e.g., MNI152NLin6Asym) ]
  [ -a | --anat_sess input : session ID of the best anatomical t1 to use (default is first session 01) ]
  [ -o | --output_dir input : output directory (default is /fmripre_dir/../cadica_fmriprepfoldername) ]
  [ -b | --brain_network_template input : brain network template - where the adult template of 7 brain networks is, optional ]
  [ -n | --num_ICs input : Number of ICs you want to estimate for melodic (default is automatic estimation) ]
  [ -m | --redoMELODIC Whether or not to rerun melodic even if it has already been run in the past]
  [ -h | --help print this help section here ]

  A Wrapper to run the first part of CICADA. Really just checks that all files are there in BIDS/fMRIPrep Standard Output and then calls first CICADA Script From There.
  First script runs MELODIC and calculates factors to determine signal to noise in next steps. Per subject, session, and task.
  Output default is <bids_dir>/derivatives/cicada

  Example: `basename $0` --input_dir='/fmripreppath/' --participant_id='102' --sess_id='01' --task_id='foodpics_run-01' --anat_sess='01' --output_dir='/pathtocicada_dir/' --brain_network_template='/pathtobrainnetworktemplate/' --num_ICs='30' --redoMELODIC

EOF
  exit 1
}

args=$(getopt -a -o f:p:s:t:l:a:o:b:n:mh --long input_dir:,participant_id:,sess_id:,task_id:,template_tag:,anat_sess:,output_dir:,brain_network_template:,num_ICs:,redoMELODIC,help -- "$@")
if [[ $? -gt 0 ]]; then
  usage
fi

eval set -- ${args}
while :
do
  case $1 in
    -f | --input_dir)   input_dir=$2    ; shift 2 ;;
    -p | --participant_id)   participant_id=$2    ; shift 2 ;;
    -s | --sess_id)   sess_id=$2    ; shift 2 ;;
    -t | --task_id)   task_id=$2    ; shift 2 ;;
    -l | --template_tag)   template_tag=$2    ; shift 2 ;;
    -a | --anat_sess)   anat_sess=$2    ; shift 2 ;;
    -o | --output_dir)   output_dir=$2    ; shift 2 ;;
    -b | --brain_network_template)   brain_network_template=$2    ; shift 2 ;;
    -n | --num_ICs)   num_ICs=$2    ; shift 2 ;;
    -m | --redoMELODIC)   redoMELODIC='1'    ; shift   ;;
    -h | --help)   usage    ; shift   ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    *) >&2 echo Unsupported option: $1
       usage ;;
  esac
done

if [[ $# -eq 0 ]]; then
  usage
fi

>&2 echo "input_dir   : ${input_dir}"
>&2 echo "participant_id    : ${participant_id} "
>&2 echo "sess_id   : ${sess_od}"
>&2 echo "task_id   : ${task_id}"
>&2 echo "template_tag   : ${tamplate_tag}"
>&2 echo "anat_sess   : ${anat_sess}"
>&2 echo "output_dir   : ${output_dir}"
>&2 echo "template_dir   : ${template_dir}"
>&2 echo "num_ICs   : ${num_ICs}"
>&2 echo "redoMELODIC   : ${redoMELODIC}"
>&2 echo "Parameters remaining are: $@"
exit 0

# OK, now check for issues with some of the defaults and reset if needed


# example of extracting tr from files for later use: This does give a long form of the number (e.g., 2.0000)
#awk '/pixdim4/{print $2}' <<< $(fslhd sub-106_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-CADICA_nonagg_bold.nii.gz)

## OK, first check that all necessary options are given
if [ "${input_dir}" = "" ]; then
  echo "Missing a fmriprep directory input"
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

if [ "${tr}" = "" ]; then
  echo "Missing a Repetition Time in seconds"
  usage
  exit
fi

if [ "${template_tag}" = "" ]; then
  echo "Missing a Template Space Name"
  usage
  exit
fi

# Now, check that these options all exist/work!
if [ ! -d "${input_dir}" ]; then
  echo "fmriprep input directory not found at ${input_dir}"
  echo "Exiting..."
  exit
fi
cd ${input_dir}
input_dir=$(pwd)

if [ ! -d "${input_dir}/sub-${participant_id}" ]; then
  echo "Participant Not Found At ${input_dir}sub-${input_id}"
  echo "Exiting..."
  exit
fi

input_subj_dir="${input_dir}/sub-${participant_id}"


if [ ! -d "${input_subj_dir}/ses-${sess_id}" ]; then
  echo "Session Not Found At ${input_subj_dir}/ses-${sess_id} for participant ${participant_id}"
  echo "Exiting..."
  exit
fi

input_sess_dir="${input_subj_dir}/ses-${sess_id}"

input_anatfol="${input_subj_dir}/ses-${anat_sess}/anat"

if [ ! -d "${input_anatfol}" ]; then
  echo "Did not find anatomy folder at ${input_anatfol}"
  echo "Is anatomy folder one level up (e.g., longitudinal averaging of anatomicals were run instead of picking the best one?)"
  echo "Exiting..."
  exit
fi

# If brain network template is not provided, assume default
if [ "${brain_network_template}" == "" ]; then
  # assume template is same location as CICADA script directory
  brain_network_template="${script_dir}/templates/cortical_subcortical_functional_atlas_guzman-Velez-2022_3mm.nii.gz"
fi

if [ ! -f ${brain_network_template} ]; then
  echo "Cannot find brain network template at ${brain_network_template}"
  exit
fi

######## WORKING FROM HERE 02/20/24 2:30PM

# Now make sure that the anatomy folder has all the necessary files
GM_anatfileid="space-${temp}_res-0${res}_label-GM_probseg.nii." # something specific in naming to the GM anatfile of interest
WM_anatfileid="space-${temp}_res-0${res}_label-WM_probseg.nii." # something specific in naming to the WM anatfile of interest
CSF_anatfileid="space-${temp}_res-0${res}_label-CSF_probseg.nii." # something specific in naming to the CSF anatfile of interest
anatmaskid="space-${temp}_res-0${res}_desc-brain_mask.nii." # something specific in naming to the anatmask of interest
anatfileid="space-${temp}_res-0${res}_desc-preproc_T1w.nii." # something specific in naming to the anatfile of interest

if ls ${input_anatfol}/*${GM_anatfileid}* 1> /dev/null 2>&1; then
  if ls ${input_anatfol}/*${WM_anatfileid}* 1> /dev/null 2>&1; then
    if ls ${input_anatfol}/*${CSFanatfileid}* 1> /dev/null 2>&1; then
      if ls ${input_anatfol}/*${anatmaskid}* 1> /dev/null 2>&1; then
        if ls ${input_anatfol}/*${anatfileid}* 1> /dev/null 2>&1; then
          echo "      Found Anatomy files"
        else
          echo "      Cannot find anatfile fitting naming scheme *${anatfileid}* for subject ${participant_id}"
          exit
        fi
      else
        echo "      Cannot find anatmask fitting naming scheme *${anatmaskid}* for subject ${participant_id}"
        exit
      fi
    else
      echo "      Cannot find CSF file fitting naming scheme *${CSFanatfileid}* for subject ${participant_id}"
      exit
    fi
  else
    echo "      Cannot find WM file fitting naming scheme *${WM_anatfileid}* for subject ${participant_id}"
    exit
  fi
else
  echo "      Cannot find GM file fitting naming scheme *${GM_anatfileid}* for subject ${participant_id}"
  exit
fi

# Now we need to look for a functional folder in the session folder

if [ ! -d "${input_sess_dir}/func" ]; then
  echo "Functional folder not found at ${input_sess_dir}/func"
  echo "Exiting..."
  exit
fi

input_func_dir="${input_sess_dir}/func"

# Now we can look to see if the functional files that we need to exist (similar to looking for the anat files)
funcmaskid="${task_id}_space-${temp}_res-0${res}_desc-brain_mask.nii." # something specific in naming to the functional mask of interest
funcfileid="${task_id}_space-${temp}_res-0${res}_desc-preproc_bold.nii." # something specific in naming to the functional file of interest
confoundfileid="${task_id}_desc-confounds_timeseries.ts"
if ls ${input_func_dir}/*${funcmaskid}* 1> /dev/null 2>&1; then
  if ls ${input_func_dir}/*${funcfileid}* 1> /dev/null 2>&1; then
    if ls ${input_func_dir}/*${confoundfileid}* 1> /dev/null 2>&1; then
      echo "      Found Functional Files"
    else
      echo "      Cannot find confound file fitting naming scheme *${confoundfileid}* for task ${task_id} for session ${sess_id} for subject ${participant_id}"
      continue
    fi
  else
    echo "      Cannot find functional file fitting naming scheme *${funcfileid}* for task ${task_id} for session ${sess_id} for subject ${participant_id}"
    continue
  fi
else
  echo "      Cannot find functional mask fitting naming scheme *${funcmaskid}* for task ${task_id} for session ${sess_id} for subject ${participant_id}"
  continue
fi

if [ "${redoMELODIC}" = "" ]; then
  redoMELODIC="0"
fi

if [ "${create_funcmask}" = "" ]; then
  create_funcmask="0"
fi

if [ "${output_dir}" = "" ]; then
  cd ${fmriprep_dir}/../
  derivatives=$(pwd)
  output_dir="${derivatives}/cadica_${temp}_res-0${res}"
fi

# let's keep logs of everything:
date
now=`date +"%Y%m%d%H%M%S"`

log_dir="${fmriprep_dir}/../logs" # directory to save log file
cadica_logfile="${log_dir}/cadica_1_output_sub_${participant_id}_session_${sess_id}_task_${task_id}_${now}.log"


(
###################### ACTUAL WORK STARTS HERE #########################################################
echo "Running for Subject ${participant_id} session ${sess_id}, task ${task_id}"

# make CADICA (Computer Assisted Denoising - ICA) folder if we have not already
if [ ! -d "${output_dir}" ]
then
  mkdir "${output_dir}"
fi
cd ${output_dir}

# make subject folder if we have not already
if [ ! -d "${output_dir}/sub-${participant_id}" ]
then
  cd ${output_dir}
  mkdir "sub-${participant_id}"
fi
output_subfol="${output_dir}/sub-${participant_id}" # main subject folder to work within
cd ${output_subfol} # go into subject folder

# make an anatmask folder on the subject level if it does not exist
if [ ! -d "${output_subfol}/anatfol" ]
then
  mkdir "${output_subfol}/anatfol"
fi
output_anatfol="${output_subfol}/anatfol"
# copy over session anatomical CSF, WM, GM masks to CADICA subject folder
find ${input_anatfol}/ -name "*${GM_anatfileid}*" -exec cp '{}' ${output_anatfol}/GM_probseg.nii.gz \;
find ${input_anatfol}/ -name "*${WM_anatfileid}*" -exec cp '{}' ${output_anatfol}/WM_probseg.nii.gz \;
find ${input_anatfol}/ -name "*${CSF_anatfileid}*" -exec cp '{}' ${output_anatfol}/CSF_probseg.nii.gz \;
find ${input_anatfol}/ -name "*${anatmaskid}*" -exec cp '{}' ${output_anatfol}/anatmask.nii.gz \;
find ${input_anatfol}/ -name "*${anatfileid}*" -exec cp '{}' ${output_anatfol}/anatfile.nii.gz \;

# point to anatomical aspects from subject space
GMprob="${output_anatfol}/GM_probseg.nii.gz"
WMprob="${output_anatfol}/WM_probseg.nii.gz"
CSFprob="${output_anatfol}/CSF_probseg.nii.gz"
anatmask="${output_anatfol}/anatmask.nii.gz"
anatfile="${output_anatfol}/anatfile.nii.gz"

echo "  Running for session ${sess_id}"

# OK, now make session folder in CADICA if it does not exist
if [ ! -d "${output_subfol}/ses-${sess_id}" ]
then
  mkdir "${output_subfol}/ses-${sess_id}"
fi

output_sess_dir="${output_subfol}/ses-${sess_id}"
cd ${output_sess_dir}

# make task folder in CADICA if we have not already
if [ ! -d "${output_sess_dir}/${task_id}" ]
then
  mkdir "${output_sess_dir}/${task_id}" # make a task folder
fi
output_task_dir="${output_sess_dir}/${task_id}"
cd ${output_task_dir}

# remake task-specific regionmask folder, just in case
if [ -d "${output_task_dir}/region_masks" ]
then
  rm -rf "region_masks"
fi
mkdir "${output_task_dir}/region_masks" # make a mask folder to put all anatomical masks of interest into
output_regionmask_dir="${output_task_dir}/region_masks"

# copy over session functional files and confounds file
find ${input_func_dir}/ -name "*${funcfileid}*" -exec cp '{}' ${output_task_dir}/funcfile_unmasked.nii.gz \;
find ${input_func_dir}/ -name "*${confoundfileid}*" -exec cp '{}' ${output_task_dir}/confounds_timeseries.csv \;

# for funcmask, use fmriprep funcmask by default, but if specified, make a new more generous mask (sometimes fMRIPrep cuts off too much)
if [ "${create_funcmask}" -eq "0" ]
then
  find ${input_func_dir}/ -name "*${funcmaskid}*" -exec cp '{}' ${output_task_dir}/funcmask.nii.gz \;
else
  echo "    Creating New Funcmask Based off fMRIPrep!"
  find ${input_func_dir}/ -name "*${funcmaskid}*" -exec cp '{}' ${output_task_dir}/funcmask_orig.nii.gz \;

fi

funcfilename="funcfile"
funcmask="${output_task_dir}/funcmask_orig.nii.gz" # func brain mask, but will relabel to get final one in next few lines
funcfile_nomask="${output_task_dir}/${funcfilename}_unmasked.nii.gz" # functional file unmasked version

# First, calculate a good resampled anatomy mask if you did not already
flirt -ref "${funcmask}" -in "${anatmask}" -out "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -usesqform -applyxfm
fslmaths "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -bin "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz"

# now find max of orig funcmask and anatmask_resam - helps make sure we do not leave out part of brain - this could be a large func mask
fslmaths "${funcmask}" -max "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" "${output_task_dir}/max_anatmask_tmp_funcmask.nii.gz"

# second, now create a smoothed anatmask and a lightly thresholded unmasked functional file, helps keep us within reasonable boundaries
fslmaths "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -s 6 -thr 0.1 -bin "${output_task_dir}/anatmask_resam_tmp_smoothed.nii.gz"
fslmaths "${funcfile_nomask}" -Tmean -thrP 5 -bin "${output_task_dir}/funcfile_unmasked_tmp_thrP5.nii.gz"

# finally, multiply the smoothed anatmask and lightly thresholded unmasked funcfile to the combined mask to get final funcmask
fslmaths "${output_task_dir}/max_anatmask_tmp_funcmask.nii.gz" -mul "${output_task_dir}/anatmask_resam_tmp_smoothed.nii.gz" -mul "${output_task_dir}/funcfile_unmasked_tmp_thrP5.nii.gz" -bin -fillh "${output_task_dir}/funcmask.nii.gz"
funcmask="${output_task_dir}/funcmask.nii.gz" # update main funcmask

# now mask your funcfile with this new funcmask
fslmaths "${funcfile_nomask}" -mul "${funcmask}" "${output_task_dir}/${funcfilename}.nii.gz"

funcfile="${output_task_dir}/${funcfilename}.nii.gz" # masked explicitely

##################### CREATE ROI MASKS ################################################################
echo "      Calculating Masks in Functional Space"

### GM, WM, and CSF Masks:
# Note: CSF mask really also includes spaces that would commonly have vessels as well.
flirt -ref "${funcmask}" -in "${GMprob}" -out "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -usesqform -applyxfm
flirt -ref "${funcmask}" -in "${WMprob}" -out "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -usesqform -applyxfm
flirt -ref "${funcmask}" -in "${CSFprob}" -out "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" -usesqform -applyxfm

# calculate susceptibility mask first then edge
fslmaths "${funcfile}" -Tmean -uthrp 20 -thr 0 -bin -mul "${funcmask}" -mul "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -bin "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz"

# Now calculate edge after removing susceptibility. Perimeter of funcmask, but not including susceptibility areas (i.e., edge must have strong enough signal)
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -ero -fmean "${output_regionmask_dir}/eroded_smoothed_tmp_func.nii.gz"
fslmaths "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -ero -fmean "${output_regionmask_dir}/eroded_smoothed_tmp_anat.nii.gz"
fslmaths ${funcmask} -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -sub "${output_regionmask_dir}/eroded_smoothed_tmp_func.nii.gz" -thr 0 "${output_regionmask_dir}/Edge_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/Edge_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/Edge_prop_final.nii.gz"

### Calculate GM, WM, and CSF by subtracting out Edge (perimeter) and susceptibility
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prop_final.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/GM_prop.nii.gz"
fslmaths "${output_regionmask_dir}/GM_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/GM_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prop_final.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/WM_prop.nii.gz"
fslmaths "${output_regionmask_dir}/WM_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/WM_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prop_final.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/CSF_prop.nii.gz"
fslmaths "${output_regionmask_dir}/CSF_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/CSF_prop_final.nii.gz"

# Can calculate an "inner CSF" by multiplying by an eroded anatmask_resam -- this may be helpful for a more target inner CSF measure. This could be in contrast to OutbrainOnly
fslmaths "${output_regionmask_dir}/anatmask_tmp_resam.nii.gz" -ero -mul ${funcmask} -bin "${output_regionmask_dir}/anatmask_resam_tmp_eroded.nii.gz"
fslmaths "${output_regionmask_dir}/anatmask_resam_tmp_eroded.nii.gz" -mul "${output_regionmask_dir}/CSF_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/InnerCSF_prop_final.nii.gz"

# Calculate WMorCSF, GMorCSF, and GMorWM regions
fslmaths "${output_regionmask_dir}/WM_tmp_prop.nii.gz" -add "${output_regionmask_dir}/CSF_tmp_prop.nii.gz" -thr 0 "${output_regionmask_dir}/WMorCSF_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/WMorCSF_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/WMorCSF_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/GM_tmp_prop.nii.gz" -add "${output_regionmask_dir}/CSF_tmp_prop.nii.gz" -thr 0 "${output_regionmask_dir}/GMorCSF_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/GMorCSF_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/GMorCSF_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/GM_tmp_prop.nii.gz" -add "${output_regionmask_dir}/WM_tmp_prop.nii.gz" -thr 0 "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/GMorWM_prop_final.nii.gz"

# The aboves are or statements, but we need "and" statements for specific boundaries
fslmaths "${output_regionmask_dir}/WMorCSF_prop_final.nii.gz" -sub "${output_regionmask_dir}/WM_prop_final.nii.gz" -sub "${output_regionmask_dir}/CSF_prop_final.nii.gz" -thr 0 "${output_regionmask_dir}/WMandCSF_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/GMorCSF_prop_final.nii.gz" -sub "${output_regionmask_dir}/GM_prop_final.nii.gz" -sub "${output_regionmask_dir}/CSF_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/GMandCSF_prop_final.nii.gz"
fslmaths "${output_regionmask_dir}/GMorWM_prop_final.nii.gz" -sub "${output_regionmask_dir}/GM_prop_final.nii.gz" -sub "${output_regionmask_dir}/WM_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/GMandWM_prop_final.nii.gz"

# Subependymal is WMandCSF plus a dilation into white matter
# first, dilate into WM
fslmaths "${output_regionmask_dir}/WMandCSF_prop_final.nii.gz" -dilM -mul "${output_regionmask_dir}/WM_prop_final.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/WMandCSF_tmp_WMdilation.nii.gz"
# now combine them:
fslmaths "${output_regionmask_dir}/WMandCSF_tmp_WMdilation.nii.gz" -add "${output_regionmask_dir}/WMandCSF_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/Subepe_prop_final.nii.gz"

# we can make WM final mask more accurate for signal now by removing Subepe from it
fslmaths "${output_regionmask_dir}/WM_prop_final.nii.gz" -sub "${output_regionmask_dir}/Subepe_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/WM_prop_final.nii.gz"

# Because WM around GM may have true BOLD signal, we can be more generous with GM by including GMWM overlap - this might be most indicative of signal with low chance of noise
fslmaths "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz" -sub "${output_regionmask_dir}/WM_prop_final.nii.gz" -sub "${output_regionmask_dir}/Subepe_prop_final.nii.gz" -mul ${funcmask} -thr 0 "${output_regionmask_dir}/GMWMlenient_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/GMWMlenient_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/GMWMlenient_prop_final.nii.gz"

# also make versions that do not worry about removing edge and susceptibility, in case you want those later
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" "${output_regionmask_dir}/Anatprob_tmp_resam.nii.gz"
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" "${output_regionmask_dir}/Inbrainprob_tmp_resam.nii.gz"

### Calculate a Not GM mask (like global signal, but without GM)
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/NotGM_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/NotGM_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/NotGM_prop_final.nii.gz"

# Also calculate a not GM OR WM mask
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/NotGMorWM_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/NotGMorWM_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/NotGMorWM_prop_final.nii.gz"

### Calculate Inbrain & Outbrain: funcmask outside of anatomical inbrain+CSF+WMCSF. Inbrain should not include WMCSF or susceptibility, since these regions cannot be trusted.
# inbrain is simply GMorWM (minus Subepe), especially since this already removes edge and susceptibility and anatomy outside of funcmask
fslmaths "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz" -thr 0.67 -bin -sub "${output_regionmask_dir}/Subepe_prop_final.nii.gz" -thr 0 -bin "${output_regionmask_dir}/Inbrain_prop_final.nii.gz"
# outbrain includes, but is not limited to, CSF. Need to remove GMorWM from funcmask. Outbrain only will also remove edge and susceptibility
# all of outbrain will include edge, susceptibility, and subepe
fslmaths "${output_regionmask_dir}/GMorWM_tmp_prop.nii.gz" -thr 0.33 "${output_regionmask_dir}/InbrainMainly_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/InbrainMainly_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/InbrainMainly_prop_final.nii.gz"

fslmaths "${funcmask}" -sub "${output_regionmask_dir}/InbrainMainly_tmp_prop.nii.gz" -thr 0 "${output_regionmask_dir}/Outbrain_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/Outbrain_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/Outbrain_prop_final.nii.gz"

# outbrainonly ideally does not include susceptibility, edge, csf, or Subepe. We can approximate this pretty well!
# This makes it so that full funcmask coverage should include: Edge + Suscept + OutbrainOnly + CSF + Subepe + Inbrain
fslmaths "${output_regionmask_dir}/Outbrain_prop_final.nii.gz" -sub "${output_regionmask_dir}/CSF_prop_final.nii.gz" -sub "${output_regionmask_dir}/Subepe_prop_final.nii.gz" -sub "${output_regionmask_dir}/Edge_prop_final.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prop_final.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/OutbrainOnly_tmp_prop.nii.gz"
fslmaths "${output_regionmask_dir}/OutbrainOnly_tmp_prop.nii.gz" -thr 0.67 -bin "${output_regionmask_dir}/OutbrainOnly_prop_final.nii.gz"


echo "      Functional Space Masks are Computed"
# relabel these for easier future coding
Edgemask="${output_regionmask_dir}/Edge_prop_final.nii.gz"
GMmask="${output_regionmask_dir}/GM_prop_final.nii.gz"
WMmask="${output_regionmask_dir}/WM_prop_final.nii.gz"
CSFmask="${output_regionmask_dir}/CSF_prop_final.nii.gz"
InnerCSFmask="${output_regionmask_dir}/InnerCSF_prop_final.nii.gz"
Outbrainmask="${output_regionmask_dir}/Outbrain_prop_final.nii.gz"
OutbrainOnlymask="${output_regionmask_dir}/OutbrainOnly_prop_final.nii.gz"
Susceptmask="${output_regionmask_dir}/Susceptibility_prop_final.nii.gz"
Inbrainmask="${output_regionmask_dir}/Inbrain_prop_final.nii.gz"
Signalmask="${output_regionmask_dir}/GMWMlenient_prop_final.nii.gz"
Subepemask="${output_regionmask_dir}/Subepe_prop_final.nii.gz"
WMCSFmask="${output_regionmask_dir}/WMandCSF_prop_final.nii.gz"
GMCSFmask="${output_regionmask_dir}/GMandCSF_prop_final.nii.gz"
GMWMmask="${output_regionmask_dir}/GMandWM_prop_final.nii.gz"
NotGMmask="${output_regionmask_dir}/NotGM_prop_final.nii.gz"

echo "      Copying of Files & Calculating Masks is Done! Now Calculating Relevant TimeSeries!"
##############################################################################################################
fslmeants -i ${funcfile} -o NotGM_${funcfilename}_timeseries.txt -m ${NotGMmask}
fslmeants -i ${funcfile} -o WM_${funcfilename}_timeseries.txt -m ${WMmask}
fslmeants -i ${funcfile} -o CSF_${funcfilename}_timeseries.txt -m ${CSFmask}
fslmeants -i ${funcfile} -o InnerCSF_${funcfilename}_timeseries.txt -m ${InnerCSFmask}
fslmeants -i ${funcfile} -o Suscept_${funcfilename}_timeseries.txt -m ${Susceptmask}
fslmeants -i ${funcfile} -o Outbrain_${funcfilename}_timeseries.txt -m ${Outbrainmask}
fslmeants -i ${funcfile} -o OutbrainOnly_${funcfilename}_timeseries.txt -m ${OutbrainOnlymask}
fslmeants -i ${funcfile} -o Edge_${funcfilename}_timeseries.txt -m ${Edgemask}
fslmeants -i ${funcfile} -o WMCSF_${funcfilename}_timeseries.txt -m ${WMCSFmask}
fslmeants -i ${funcfile} -o Subepe_${funcfilename}_timeseries.txt -m ${Subepemask}

echo "      Relevant TimeSeries Are Computed! Now Running MELODIC!"
##################### MELODIC ICA DECOMPOSITION SECTION ######################################################
# let's finally get to melodic IC decomposition
cd ${output_task_dir}

# MELODIC SECTION
# make session folder if we have not already
melfol="${output_task_dir}/melodic"
if [ -d "${output_task_dir}/melodic" ]
then
  if [ "${redoMELODIC}" -eq "0" ]
  then
    echo "    Not rerunning Melodic"
  else
    cd ${output_task_dir}
    rm -rf melodic
    mkdir melodic

    if [ "${num_ICs}" = "" ]; then
      echo "      Running: melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
      melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
    else
      echo "      Running: melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --dim="${num_ICs}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
      melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --dim="${num_ICs}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
    fi


  fi
fi

if [ ! -d "${output_task_dir}/melodic" ]
then
  mkdir "${output_task_dir}/melodic"
  if [ "${num_ICs}" = "" ]; then
    echo "      Running: melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
    melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
  else
    echo "      Running: melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --dim="${num_ICs}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
    melodic --in="${funcfile}" --outdir=${melfol} --mask="${funcmask}" --dim="${num_ICs}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
  fi
fi


# Merge the probability maps and threshold maps which may prove useful later
probmapnames="$(ls ${melfol}/stats/probmap_* | sort -V)"
fslmerge -t ${melfol}/ICprobabilities.nii.gz ${probmapnames}

# Calculate a 99% ICprobabilities mask
fslmaths "${melfol}/ICprobabilities.nii.gz" -thr 0.99 -bin "${melfol}/ICprobabilities_99percent.nii.gz"

thresholdnames="$(ls ${melfol}/stats/thresh_zstat* | sort -V)"
fslmerge -t ${melfol}/ICthresh_zstat.nii.gz ${thresholdnames}

cd ${output_task_dir}

echo "      Melodic is Complete! Now Calculating Smoothness and Probabilities"
#########################################################################################
## Calculate smoothness retention before and after 6mm gaussian smoothing of IC
# Get explained variance
cd ${output_task_dir}
if [ -d "ROIcalcs" ]
then
  rm -rf ROIcalcs
fi
mkdir "ROIcalcs"
ROIcalcfol="${output_task_dir}/ROIcalcs"

# grab explained variance percentage in case of use in weighting the Data
awk '{print $1}' ${melfol}/melodic_ICstats > ${ROIcalcfol}/IC_exp_variance.txt

# Calculate the higher probabilities to use
fslmaths ${melfol}/ICprobabilities.nii.gz -thr 0.95 "${ROIcalcfol}/highprob_tmp_prop.nii.gz"
fslmaths "${ROIcalcfol}/highprob_prop.nii.gz" -bin "${ROIcalcfol}/highprob_tmp_mask.nii.gz"

# calculate a nonthresholded version of the z stats, absolute valued
fslmaths ${melfol}/ICthresh_zstat.nii.gz -abs "${ROIcalcfol}/fullvolICA_tmp_nothresh.nii.gz"

# calculate a smoothed 6mm gauss version of nonthresholded too
fslmaths ${melfol}/ICthresh_zstat.nii.gz -s 6 -abs "${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz"

# Calculate what you need form these to compare before and after smoothing to get a smoothing retention
fslstats -t "${ROIcalcfol}/fullvolICA_tmp_nothresh.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_nothresh_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_nothresh_ICnumvoxels.txt

fslstats -t "${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_smoothed_nothresh_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_smoothed_nothresh_ICnumvoxels.txt

#########################################################################################
# Now set up brain network template for easy comparisons
echo "      Smoothing Checks are Complete! Now Setting Up Brain Network Template"
# check if the resampled template already exists, if not, then make it
template_dir="${output_dir}/templates"
if [ ! -d "${template_dir}" ]; then
  mkdir "${template_dir}"
fi

if [ ! -f "${template_dir}/network_template_${task_id}.nii.gz" ]; then
  echo "        Template Resampled Not Found, Resampling one now at ${template_dir}/network_template_${task_id}.nii.gz."
  flirt -ref "${funcmask}" -in "${brain_network_template}" -out "${template_dir}/network_template_${task_id}.nii.gz" -usesqform -applyxfm -interp nearestneighbour
fi
# OK, now we have a resampled template that we can use, should have 1-7 values
# 1: Medial Visual, 2: Sensory Motor, 3: Dorsal Attention, 4: Ventral Attention, 5: FrontoParietal, 6: Default Mode Network, 7: Subcortical
network_template="${template_dir}/network_template_${task_id}.nii.gz"
#########################################################################################
# NOW, do all the relevant calculations!
echo "      Brain Network Template Calculated. Now finally do all the relevant calculations!"
calcfile="${ROIcalcfol}/highprob_tmp_prop.nii.gz"

# GM voxels
fslmaths "${calcfile}" -mul "${GMmask}" "${ROIcalcfol}/GMICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/GMICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/GM_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/GM_ICnumvoxels.txt

# GM WM lenient voxels, labeled as "signal voxels"
fslmaths "${calcfile}" -mul "${Signalmask}" "${ROIcalcfol}/SignalICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/SignalICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Signal_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Signal_ICnumvoxels.txt

# WM voxels
fslmaths "${calcfile}" -mul "${WMmask}" "${ROIcalcfol}/WMICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/WMICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/WM_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/WM_ICnumvoxels.txt

# CSF voxels
fslmaths "${calcfile}" -mul "${CSFmask}" "${ROIcalcfol}/CSFICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/CSFICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/CSF_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/CSF_ICnumvoxels.txt

# Inner CSF voxels
fslmaths "${calcfile}" -mul "${InnerCSFmask}" "${ROIcalcfol}/InnerCSFICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/InnerCSFICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/InnerCSF_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/InnerCSF_ICnumvoxels.txt

# Edge voxels
fslmaths "${calcfile}" -mul "${Edgemask}" "${ROIcalcfol}/EdgeICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/EdgeICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Edge_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Edge_ICnumvoxels.txt

# Outbrain voxels including CSF
fslmaths "${calcfile}" -mul "${Outbrainmask}" "${ROIcalcfol}/OutbrainICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/OutbrainICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Outbrain_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Outbrain_ICnumvoxels.txt

# Outbrain only voxels (no inner CSF, edge, or suscept) - good for sinuses
fslmaths "${calcfile}" -mul "${OutbrainOnlymask}" "${ROIcalcfol}/OutbrainOnlyICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/OutbrainOnlyICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/OutbrainOnly_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/OutbrainOnly_ICnumvoxels.txt

# Inbrain voxels
fslmaths "${calcfile}" -mul "${Inbrainmask}" "${ROIcalcfol}/InbrainICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/InbrainICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Inbrain_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Inbrain_ICnumvoxels.txt

# WM CSF Boundary + surrounding WM voxels (Subependymal)
fslmaths "${calcfile}" -mul "${Subepemask}" "${ROIcalcfol}/SubepeICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/SubepeICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Subepe_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Subepe_ICnumvoxels.txt

# Susceptibility voxels
fslmaths "${calcfile}" -mul "${Susceptmask}" "${ROIcalcfol}/SusceptICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/SusceptICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Suscept_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/Suscept_ICnumvoxels.txt

# WMCSF Boundary voxels
fslmaths "${calcfile}" -mul "${WMCSFmask}" "${ROIcalcfol}/WMCSFICA_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/WMCSFICA_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/WMCSF_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/WMCSF_ICnumvoxels.txt

# Now also grab information for the 7 networks:
# 1: Medial Visual, 2: Sensory Motor, 3: Dorsal Attention, 4: Ventral Attention, 5: FrontoParietal, 6: Default Mode Network, 7: Subcortical
tag="MedialVisual"
val="1"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="SensoryMotor"
val="2"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="DorsalAttention"
val="3"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="VentralAttention"
val="4"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="FrontoParietal"
val="5"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="DefaultModeNetwork"
val="6"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

tag="Subcortical"
val="7"
fslmaths "${network_template}" -uthr ${val} -thr ${val} -bin -mul "${Signalmask}" -mul "${calcfile}" "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz"
fslstats -t "${ROIcalcfol}/${tag}_tmp_weighted.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICmean.txt"
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > "${ROIcalcfol}/${tag}_ICnumvoxels.txt"

# remove extraneous files that will no longer be needed for future work (reduce size of total output)
rm ${ROIcalcfol}/curr_tmp_calc.txt
rm -f ${output_regionmask_dir}/*_tmp_*.nii.gz
rm -f ${output_task_dir}/*_tmp_*.nii.gz
rm -f ${ROIcalcfol}/*_tmp_*.nii.gz
###########################################################################################################

echo "    Done with task ${task_id} for session ${sess_id} for subject ${participant_id}"
) 2>&1 | tee ${cadica_logfile}
