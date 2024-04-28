#!/bin/bash

# Run first script/function of CICADA
# NOTE: make sure your version of FSl is creating .nii.gz files before opening Matlab and running any scripts.
# e.g.,: in bash_profile:
# FSLOUTPUTTYPE=NIFTI_GZ
# export FSLOUTPUTTYPE


# immediately exit upon common error
set -euo pipefail
echo

usage(){
>&2 cat << EOF
Usage: `basename $0` -o /cicada_dir/sub_id/ses_id/task_id -F funcfile.nii.gz -f funcmask.nii.gz -C confoundsfile.csv -A anatfile.nii.gz -a anatmask.nii.gz -g gm_prob.nii.gz -w wm_prob.nii.gz -c csf_prob.nii.gz
  [ -o input : (necessary) output_dir: where to put the outputs, should be cicada_dir/sub_id/sess_id/task_id. Make it two levels deeper than subject_id folder! ]
  [ -F input : (necessary) funcfile: unmasked functional file ]
  [ -f input : (necessary) funcmask: current functional mask ]
  [ -C input : (necessary) confoundsfile: csv of confounds as columns with headers. Need 6 motion parameters (e.g., rot_x, trans_x, etc.), dvars, framewise_displacement, csf, white_matter, and global_signal ]
  [ -A input : (preferred) anatfile: anatomy file, will default to MNI 2009c asym T1 ]
  [ -a input : (preferred) anatmask: anatomy mask, will default to MNI 2009c asym mask ]
  [ -g input : (preferred) gm_prob: grey matter probability file, will default to MNI 2009c asym ]
  [ -w input : (preferred) wm_prob: white matter probability file, will default to MNI 2009c asym ]
  [ -c input : (preferred) csf_prob: CSF probability file, will default to MNI ]
  [ -m input : (optional) mel_fol: melodic folder to use (instead of redoing melodic) ]
  [ -h : help: print this help section here ]

  Run the first part of CICADA (make masks, run MELODIC, and do relevant calculations to prepare for IC selection).
  Important Notes:
  1. Make sure a subject folder is 2 levels above the output_dir. Strongly recommended that output dir be cicada_dir/subj_dir/sess_dir/task_dir
  2. Confounds file must be 1 confound per column, comma separated, with first row as a header/labels for each confound. fMRIPrep's confounds.tsv file works well for this.

EOF
  exit 1
}

if [ $# -lt 8 -o $# -gt 21 ]; then
  echo "Number of inputs ($#) are incorrect"
  echo
	usage
  echo
  echo
fi

# initialize things
# Initialize things if not specified ("x" will say "does not exist")
output_dir="x"
funcfile="x"
funcmask="x"
confoundsfile="x"
anatfile="x"
anatmask="x"
csf_prob="x"
wm_prob="x"
gm_prob="x"
brain_network_template="x"
mel_fol="x"
cd "$(dirname $0)"
script_dir="$(pwd)"
cd "${script_dir}/../templates"
template_dir="$(pwd)"
mni_dir="${template_dir}/mni_icbm152_nlin_asym_09c"
cd "${script_dir}"

while getopts "o:F:f:C:A:a:g:w:c:m:h" opt; do
    case "${opt}" in
      o)
    		output_dir=${OPTARG}
    		;;
      F)
        funcfile=${OPTARG}
        ;;
      f)
        funcmask=${OPTARG}
        ;;
      C)
        confoundsfile=${OPTARG}
        ;;
      A)
        anatfile=${OPTARG}
        ;;
      a)
        anatmask=${OPTARG}
        ;;
      g)
        gm_prob=${OPTARG}
        ;;
      w)
        wm_prob=${OPTARG}
        ;;
      c)
        csf_prob=${OPTARG}
        ;;
      m)
        mel_fol=${OPTARG}
        ;;
      h)
        usage
        ;;
      *)
        usage
        ;;
    esac
done

# Now, if certain inputs are not given, either through an error or make them into their defaults if they have them:
if [ "${output_dir}" = "x" ]; then
  echo "No output directory specified"
  usage
  exit
fi
if [ "${funcfile}" = "x" ]; then
  echo "Missing an unmasked functional file"
  usage
  exit
fi
if [ "${funcmask}" = "x" ]; then
  echo "Missing a functional mask file"
  usage
  exit
fi
if [ "${confoundsfile}" = "x" ]; then
  echo "Missing a confounds file. Should be .csv, confounds are in columns with headers."
  usage
  exit
fi
if [ "${anatfile}" = "x" ]; then
  echo "Using MNI anatfile"
  anatfile="${mni_dir}/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz"
fi
if [ "${anatmask}" = "x" ]; then
  echo "Using MNI anatmask"
  anatmask="${mni_dir}/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz"
fi
if [ "${gm_prob}" = "x" ]; then
  echo "Using MNI gm_prob"
  gm_prob="${mni_dir}/mni_icbm152_gm_tal_nlin_asym_09c.nii.gz"
fi
if [ "${wm_prob}" = "x" ]; then
  echo "Using MNI wm_prob"
  wm_prob="${mni_dir}/mni_icbm152_wm_tal_nlin_asym_09c.nii.gz"
fi
if [ "${csf_prob}" = "x" ]; then
  echo "Using MNI csf_prob"
  csf_prob="${mni_dir}/mni_icbm152_csf_tal_nlin_asym_09c.nii.gz"
fi
if [ "${brain_network_template}" == "x" ]; then
  # assume template is same location as CICADA script directory
  brain_network_template="${template_dir}/cortical_subcortical_functional_atlas_guzman-Velez-2022_3mm.nii.gz"
fi
# mel_fol being set to x tells it to re run it regardless
if [ "${mel_fol}" = "x" ]; then
  mel_fol="${output_dir}/melodic"
  echo "Melodic folder not found. Will be run and stored in output directory: ${mel_fol}"
  # if this directory already existed, remove it
  if [ -d "${mel_fol}" ]; then
    echo "Remnants of a Melodic folder found. Removing it."
    rm -rf "${mel_fol}"
  fi
fi
echo
# Finally, list what you have
>&2 echo "output_dir: ${output_dir}"
>&2 echo "funcfile: ${funcfile} "
>&2 echo "funcmask: ${funcmask}"
>&2 echo "confoundsfile: ${confoundsfile}"
>&2 echo "anatfile: ${anatfile}"
>&2 echo "anatmask: ${anatmask}"
>&2 echo "gm_prob: ${gm_prob}"
>&2 echo "wm_prob: ${wm_prob}"
>&2 echo "csf_prob: ${csf_prob}"
>&2 echo "mel_fol: ${mel_fol}"
>&2 echo "script_dir: ${script_dir}"
echo

# Now check that these options all exist/work!
if [ ! -d "${template_dir}" ]; then
  echo "Cannot find template directory at ${template_dir}"
  exit
fi

if [ ! -d "${mni_dir}" ]; then
  echo "Cannot find mni directory at ${mni_dir}"
  exit
fi

if [ ! -f ${funcfile} ]; then
  echo "Cannot find functional file at ${funcfile}"
  exit
fi

if [ ! -f ${funcmask} ]; then
  echo "Cannot find functional mask at ${funcmask}"
  exit
fi

if [ ! -f ${confoundsfile} ]; then
  echo "Cannot find confoundsfile at ${confoundsfile}"
  exit
fi

if [ ! -f ${anatfile} ]; then
  echo "Cannot find anatfile at ${anatfile}"
  exit
fi

if [ ! -f ${anatmask} ]; then
  echo "Cannot find anatmask at ${anatmask}"
  exit
fi

if [ ! -f ${gm_prob} ]; then
  echo "Cannot find grey matter probability file at ${gm_prob}"
  exit
fi

if [ ! -f ${wm_prob} ]; then
  echo "Cannot find white matter probability file at ${wm_prob}"
  exit
fi

if [ ! -f ${csf_prob} ]; then
  echo "Cannot find CSF probability file at ${csf_prob}"
  exit
fi

if [ ! -f ${brain_network_template} ]; then
  echo "Cannot find brain network template at ${brain_network_template}"
  exit
fi

# if output dir does not exist, make it (do not delete the whole thing because it may contain a melodic folder or log file that you want)
if [ ! -d "${output_dir}" ]; then
  mkdir -p ${output_dir}
fi


# let's keep logs of everything in this script:
echo
echo "Running First CICADA Script for ${output_dir}!"
date
now=`date +"%Y%m%d%H%M%S"`

log_dir="${output_dir}/logs" # directory to save log file
if [ ! -d "${log_dir}" ]; then
  mkdir ${log_dir}
fi
cicada_logfile="${log_dir}/cicada_1.log"


(
###################### ACTUAL WORK STARTS HERE #########################################################

# ses_dir is one level up, Subject folder should be two levels up, cicada_dir should be three
cd ${output_dir}/../
ses_dir="$(pwd)"
cd ../
subj_dir="$(pwd)"
cd ../
cicada_home_dir="$(pwd)"
cd ${output_dir}
task_dir="$(pwd)"
# get just the task name
task_name="$(basename ${PWD})"


# Make an anatomy directory for reference, one levels up (should be session folder)
cd ${output_dir}

# make an anatmask folder on the subject level if it does not exist
if [ ! -d "${ses_dir}/anatfol" ]
then
  mkdir "${ses_dir}/anatfol"
fi
output_anatfol="${ses_dir}/anatfol"

# Copy over files for easy future reference
cp "${gm_prob}" "${output_anatfol}/GM_probseg.nii.gz"
cp "${wm_prob}" "${output_anatfol}/WM_probseg.nii.gz"
cp "${csf_prob}" "${output_anatfol}/CSF_probseg.nii.gz"
cp "${anatfile}" "${output_anatfol}/anatfile.nii.gz"
cp "${anatmask}" "${output_anatfol}/anatmask.nii.gz"

# point to anatomical aspects from subject space
GMprob="${output_anatfol}/GM_probseg.nii.gz"
WMprob="${output_anatfol}/WM_probseg.nii.gz"
CSFprob="${output_anatfol}/CSF_probseg.nii.gz"
anatmask="${output_anatfol}/anatmask.nii.gz"
anatfile="${output_anatfol}/anatfile.nii.gz"


# copy over session functional files and confounds file
funcfilename="funcfile"
cp "${funcfile}" "${output_dir}/${funcfilename}_unmasked.nii.gz"
cp "${funcmask}" "${output_dir}/funcmask_orig.nii.gz"
cp "${confoundsfile}" "${output_dir}/confounds_timeseries.csv"

# remake task-specific regionmask folder, just in case
if [ -d "${output_dir}/region_masks" ]
then
  rm -rf "region_masks"
fi
mkdir "${output_dir}/region_masks" # make a mask folder to put all anatomical masks of interest into
output_regionmask_dir="${output_dir}/region_masks"

echo "    Creating New Funcmask To Include Outer Areas Of Brain!"
funcmask="${output_dir}/funcmask_orig.nii.gz" # func brain mask, but will relabel to get final one in next few lines
funcfile_nomask="${output_dir}/${funcfilename}_unmasked.nii.gz" # functional file unmasked version

# create an eroded anatmask which could be helpful later for looking more towards the center of the brain
# four erodes is likely good for standard 2mm space
fslmaths "${anatmask}" -ero -ero -ero -ero -ero "${output_anatfol}/anatmask_eroded.nii.gz"
anatmask_eroded="${output_anatfol}/anatmask_eroded.nii.gz"

# First, calculate a good resampled anatomy mask if you did not already. Do not make anatmask_resam temporary because potentially useful later.
flirt -ref "${funcmask}" -in "${anatmask}" -out "${output_regionmask_dir}/anatmask_resam.nii.gz" -usesqform -applyxfm
fslmaths "${output_regionmask_dir}/anatmask_resam.nii.gz" -bin "${output_regionmask_dir}/anatmask_resam.nii.gz"

# Do it for the eroded one too!
flirt -ref "${funcmask}" -in "${anatmask_eroded}" -out "${output_regionmask_dir}/anatmask_eroded_resam.nii.gz" -usesqform -applyxfm
fslmaths "${output_regionmask_dir}/anatmask_eroded_resam.nii.gz" -bin "${output_regionmask_dir}/anatmask_eroded_resam.nii.gz"

# now find max of orig funcmask and anatmask_resam - helps make sure we do not leave out part of brain - this could be a large func mask
fslmaths "${funcmask}" -max "${output_regionmask_dir}/anatmask_resam.nii.gz" "${output_dir}/max_anatmask_tmp_funcmask.nii.gz"

# second, now create a smoothed anatmask and a lightly thresholded unmasked functional file, helps keep us within reasonable boundaries
fslmaths "${output_regionmask_dir}/anatmask_resam.nii.gz" -s 6 -thr 0.1 -bin "${output_dir}/anatmask_resam_tmp_smoothed.nii.gz"
fslmaths "${funcfile_nomask}" -Tmean -thrP 5 -bin "${output_dir}/funcfile_unmasked_tmp_thrP5.nii.gz"

# finally, multiply the smoothed anatmask and lightly thresholded unmasked funcfile to the combined mask to get final funcmask
fslmaths "${output_dir}/max_anatmask_tmp_funcmask.nii.gz" -mul "${output_dir}/anatmask_resam_tmp_smoothed.nii.gz" -mul "${output_dir}/funcfile_unmasked_tmp_thrP5.nii.gz" -bin -fillh "${output_dir}/funcmask.nii.gz"
funcmask="${output_dir}/funcmask.nii.gz" # update main funcmask

# now mask your unmasked funcfile with this new funcmask
fslmaths "${funcfile_nomask}" -mul "${funcmask}" "${output_dir}/${funcfilename}.nii.gz"

funcfile="${output_dir}/${funcfilename}.nii.gz" # masked explicitely

##################### CREATE ROI MASKS ################################################################
echo "    Calculating ROI Masks in Functional Space"

### GM, WM, and CSF Masks:
# Note: CSF mask really also includes spaces that would commonly have vessels as well.
flirt -ref "${funcmask}" -in "${GMprob}" -out "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -usesqform -applyxfm
flirt -ref "${funcmask}" -in "${WMprob}" -out "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -usesqform -applyxfm
flirt -ref "${funcmask}" -in "${CSFprob}" -out "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" -usesqform -applyxfm

# calculate susceptibility mask first then edge
fslmaths "${funcfile}" -Tmean -thr 0 -mul "${output_regionmask_dir}/anatmask_resam.nii.gz" -mul "${funcmask}" "${output_regionmask_dir}/funcfile_tmp_tmean_thresholded.nii.gz"
range_vals=($(fslstats ${output_regionmask_dir}/funcfile_tmp_tmean_thresholded.nii.gz -l 0.01 -r))
fslmaths "${output_regionmask_dir}/funcfile_tmp_tmean_thresholded.nii.gz" -div "${range_vals[1]}" "${output_regionmask_dir}/Suscept_almost_tmp_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Suscept_almost_tmp_prob.nii.gz" -mul -1 -add 0.5 -thr 0 -mul 2 -mul "${funcmask}" -mul "${output_regionmask_dir}/anatmask_resam.nii.gz" "${output_regionmask_dir}/Susceptibility_tmp_prob.nii.gz" # keep it within brain anatomy, and also weird math so below 0.5 becomes 0
range_vals=($(fslstats ${output_regionmask_dir}/Susceptibility_tmp_prob.nii.gz -R)) # maximum value is OK here
fslmaths "${output_regionmask_dir}/Susceptibility_tmp_prob.nii.gz" -div "${range_vals[1]}" "${output_regionmask_dir}/Susceptibility_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/Susceptibility_mask.nii.gz"
#fslmaths "${funcfile}" -Tmean -uthrp 20 -thr 0 -bin -mul "${funcmask}" -mul "${output_regionmask_dir}/anatmask_resam.nii.gz" -bin "${output_regionmask_dir}/Susceptibility_mask.nii.gz"

# Now calculate edge after removing susceptibility. Perimeter of funcmask, but not including susceptibility areas (i.e., edge must have strong enough signal)
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/Susceptibility_mask.nii.gz" -ero -fmean "${output_regionmask_dir}/eroded_smoothed_tmp_func.nii.gz"
fslmaths "${output_regionmask_dir}/anatmask_resam.nii.gz" -ero -fmean "${output_regionmask_dir}/eroded_smoothed_tmp_anat.nii.gz"
fslmaths ${funcmask} -sub "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -sub "${output_regionmask_dir}/eroded_smoothed_tmp_func.nii.gz" -thr 0 "${output_regionmask_dir}/Edge_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Edge_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/Edge_mask.nii.gz"

### Calculate GM, WM, and CSF by subtracting out Edge (perimeter) and susceptibility
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prob.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/GM_prob.nii.gz"
fslmaths "${output_regionmask_dir}/GM_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/GM_mask.nii.gz"
fslmaths "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prob.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/WM_prob.nii.gz"
fslmaths "${output_regionmask_dir}/WM_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/WM_mask.nii.gz"
fslmaths "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" -sub "${output_regionmask_dir}/Edge_prob.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -thr 0 -mul "${funcmask}" "${output_regionmask_dir}/CSF_prob.nii.gz"
fslmaths "${output_regionmask_dir}/CSF_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/CSF_mask.nii.gz"

# Make an extended GM_prob which goes more into outbrain, just so this can be more taken into account
fslmaths "${output_regionmask_dir}/GM_prob.nii.gz" -bin "${output_regionmask_dir}/GM_tmp_allmask.nii.gz"
fslmaths "${output_regionmask_dir}/GM_prob.nii.gz" -fmean -sub "${output_regionmask_dir}/GM_tmp_allmask.nii.gz" -thr 0 "${output_regionmask_dir}/GM_tmp_outerprob.nii.gz"
fslmaths "${output_regionmask_dir}/GM_tmp_outerprob.nii.gz" -add "${output_regionmask_dir}/GM_prob.nii.gz" "${output_regionmask_dir}/GM_extended_prob.nii.gz" # this is GM probability, but extended out a little bit. Useful for better QC testing (e.g., grab outbrain not potentially impacted by GM)

# Can calculate an "inner CSF" by multiplying by an eroded anatmask_resam -- this may be helpful for a more target inner CSF measure. This could be in contrast to OutbrainOnly
fslmaths "${output_regionmask_dir}/anatmask_eroded_resam.nii.gz" -mul "${output_regionmask_dir}/CSF_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/InnerCSF_mask.nii.gz"

# Also do inner WM, this will be useful to make subepe mask
fslmaths "${output_regionmask_dir}/anatmask_eroded_resam.nii.gz" -mul "${output_regionmask_dir}/WM_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/InnerWM_mask.nii.gz"

# OK, now get subepe by fmean both innercsf and inner wm, multiplying, and scaling
fslmaths "${output_regionmask_dir}/CSF_prob.nii.gz" -mul "${output_regionmask_dir}/InnerCSF_mask.nii.gz" -fmean "${output_regionmask_dir}/InnerCSF_tmp_smoothed.nii.gz"
fslmaths "${output_regionmask_dir}/WM_prob.nii.gz" -mul "${output_regionmask_dir}/InnerWM_mask.nii.gz" -fmean "${output_regionmask_dir}/InnerWM_tmp_smoothed.nii.gz"
fslmaths "${output_regionmask_dir}/InnerCSF_tmp_smoothed.nii.gz" -mul "${output_regionmask_dir}/InnerWM_tmp_smoothed.nii.gz" -mul 4 "${output_regionmask_dir}/Subepe_tmp_prob.nii.gz"
range_vals=($(fslstats ${output_regionmask_dir}/Subepe_tmp_prob.nii.gz -l 0.01 -r))
fslmaths "${output_regionmask_dir}/Subepe_tmp_prob.nii.gz" -div "${range_vals[1]}" "${output_regionmask_dir}/Subepe_almost_tmp_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Subepe_almost_tmp_prob.nii.gz" -thr 1 -bin "${output_regionmask_dir}/Subepe_almost_tmp_mask.nii.gz" # binarize anything above 1
fslmaths "${output_regionmask_dir}/Subepe_almost_tmp_prob.nii.gz" -uthr 1 -add "${output_regionmask_dir}/Subepe_almost_tmp_mask.nii.gz" "${output_regionmask_dir}/Subepe_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Subepe_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/Subepe_mask.nii.gz"

# Calculate WMorCSF, GMorCSF, and GMorWM regions, can be helpful
fslmaths "${output_regionmask_dir}/WM_prob.nii.gz" -add "${output_regionmask_dir}/CSF_prob.nii.gz" -thr 0 "${output_regionmask_dir}/WMorCSF_tmp_prob.nii.gz"
fslmaths "${output_regionmask_dir}/WMorCSF_tmp_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/WMorCSF_mask.nii.gz"
fslmaths "${output_regionmask_dir}/GM_prob.nii.gz" -add "${output_regionmask_dir}/CSF_prob.nii.gz" -thr 0 "${output_regionmask_dir}/GMorCSF_tmp_prob.nii.gz"
fslmaths "${output_regionmask_dir}/GMorCSF_tmp_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/GMorCSF_mask.nii.gz"
fslmaths "${output_regionmask_dir}/GM_prob.nii.gz" -add "${output_regionmask_dir}/WM_prob.nii.gz" -sub "${output_regionmask_dir}/Subepe_prob.nii.gz" -thr 0 "${output_regionmask_dir}/GMorWM_prob.nii.gz" # don't include subepe in there! It is different
fslmaths "${output_regionmask_dir}/GMorWM_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/GMorWM_mask.nii.gz"

# The aboves are or statements, but we need "and" statements for specific boundaries
fslmaths "${output_regionmask_dir}/WMorCSF_mask.nii.gz" -sub "${output_regionmask_dir}/WM_mask.nii.gz" -sub "${output_regionmask_dir}/CSF_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/WMandCSF_mask.nii.gz"
fslmaths "${output_regionmask_dir}/GMorCSF_mask.nii.gz" -sub "${output_regionmask_dir}/GM_mask.nii.gz" -sub "${output_regionmask_dir}/CSF_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/GMandCSF_mask.nii.gz"
fslmaths "${output_regionmask_dir}/GMorWM_mask.nii.gz" -sub "${output_regionmask_dir}/GM_mask.nii.gz" -sub "${output_regionmask_dir}/WM_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/GMandWM_mask.nii.gz"

# Subependymal is WMandCSF plus a dilation into white matter
# first, dilate into WM
#fslmaths "${output_regionmask_dir}/WMandCSF_mask.nii.gz" -dilM -mul "${output_regionmask_dir}/WM_mask.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/WMandCSF_tmp_WMdilation.nii.gz"
# now combine them:
#fslmaths "${output_regionmask_dir}/WMandCSF_tmp_WMdilation.nii.gz" -add "${output_regionmask_dir}/WMandCSF_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/Subepe_mask.nii.gz"

# we can make WM final mask more accurate for signal now by removing Subepe from it
fslmaths "${output_regionmask_dir}/WM_mask.nii.gz" -sub "${output_regionmask_dir}/Subepe_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/WM_adj_mask.nii.gz"

# Because WM around GM may have true BOLD signal, we can be more generous with GM by including GMWM overlap - this might be most indicative of signal with low chance of noise. Need to sub WM mask and subepe mask not prob here to maintain GMWM leniency
fslmaths "${output_regionmask_dir}/GMorWM_prob.nii.gz" -sub "${output_regionmask_dir}/WM_adj_mask.nii.gz" -sub "${output_regionmask_dir}/Subepe_mask.nii.gz" -thr 0 -mul ${funcmask} "${output_regionmask_dir}/GMWMlenient_prob.nii.gz"
fslmaths "${output_regionmask_dir}/GMWMlenient_prob.nii.gz" -thrP 67 -add "${output_regionmask_dir}/GM_mask.nii.gz" -bin "${output_regionmask_dir}/GMWMlenient_mask.nii.gz" # include adding in GM mask again in case small spots were removed from thresholding

# also make versions that do not worry about removing edge and susceptibility, in case you want those later
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/CSFprob_tmp_resam.nii.gz" "${output_regionmask_dir}/Anatprob_tmp_resam.nii.gz"
fslmaths "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -add "${output_regionmask_dir}/WMprob_tmp_resam.nii.gz" "${output_regionmask_dir}/Inbrainprob_tmp_resam.nii.gz"

### Calculate a Not GM mask (like global signal, but without GM)
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/GMprob_tmp_resam.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/NotGM_prob.nii.gz"
fslmaths "${output_regionmask_dir}/NotGM_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/NotGM_mask.nii.gz"

# Also calculate a not GM OR WM mask
fslmaths "${funcmask}" -sub "${output_regionmask_dir}/GMorWM_prob.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/NotGMorWM_prob.nii.gz"
fslmaths "${output_regionmask_dir}/NotGMorWM_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/NotGMorWM_mask.nii.gz"

### Calculate Inbrain & Outbrain: funcmask outside of anatomical inbrain+CSF+WMCSF. Inbrain should not include WMCSF or susceptibility, since these regions cannot be trusted.
# inbrain is simply GMorWM (minus Subepe), especially since this already removes edge and susceptibility and anatomy outside of funcmask
# subepe should have already been removed, but it doesn't hurt to double check
fslmaths "${output_regionmask_dir}/GMorWM_prob.nii.gz" -thrP 67 -bin -sub "${output_regionmask_dir}/Subepe_mask.nii.gz" -thr 0 -bin "${output_regionmask_dir}/Inbrain_mask.nii.gz"
# outbrain includes, but is not limited to, CSF. Need to remove GMorWM from funcmask. Outbrain only will also remove edge and susceptibility
# all of outbrain will include edge, susceptibility, and subepe
fslmaths "${output_regionmask_dir}/GMorWM_prob.nii.gz" -thrP 33 "${output_regionmask_dir}/InbrainMainly_tmp_prob.nii.gz"
fslmaths "${output_regionmask_dir}/InbrainMainly_tmp_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/InbrainMainly_mask.nii.gz"

fslmaths "${funcmask}" -sub "${output_regionmask_dir}/InbrainMainly_tmp_prob.nii.gz" -thr 0 "${output_regionmask_dir}/Outbrain_prob.nii.gz"
fslmaths "${output_regionmask_dir}/Outbrain_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/Outbrain_mask.nii.gz"

# outbrainonly ideally does not include susceptibility, edge, csf, or Subepe. We can approximate this pretty well!
# This makes it so that full funcmask coverage should include: Edge + Suscept + OutbrainOnly + CSF + Subepe + Inbrain
fslmaths "${output_regionmask_dir}/Outbrain_prob.nii.gz" -sub "${output_regionmask_dir}/CSF_prob.nii.gz" -sub "${output_regionmask_dir}/Subepe_prob.nii.gz" -sub "${output_regionmask_dir}/Edge_prob.nii.gz" -sub "${output_regionmask_dir}/Susceptibility_prob.nii.gz" -mul "${funcmask}" -thr 0 "${output_regionmask_dir}/OutbrainOnly_prob.nii.gz"
fslmaths "${output_regionmask_dir}/OutbrainOnly_prob.nii.gz" -thrP 67 -bin "${output_regionmask_dir}/OutbrainOnly_mask.nii.gz"


echo "      Functional Space Masks are Computed"
# relabel these for easier future coding
Edgemask="${output_regionmask_dir}/Edge_mask.nii.gz"
GMmask="${output_regionmask_dir}/GM_mask.nii.gz"
WMmask="${output_regionmask_dir}/WM_adj_mask.nii.gz" # so it does not include Subepe
CSFmask="${output_regionmask_dir}/CSF_mask.nii.gz"
InnerCSFmask="${output_regionmask_dir}/InnerCSF_mask.nii.gz"
Outbrainmask="${output_regionmask_dir}/Outbrain_mask.nii.gz"
OutbrainOnlymask="${output_regionmask_dir}/OutbrainOnly_mask.nii.gz"
Susceptmask="${output_regionmask_dir}/Susceptibility_mask.nii.gz"
Inbrainmask="${output_regionmask_dir}/Inbrain_mask.nii.gz"
Signalmask="${output_regionmask_dir}/GMWMlenient_mask.nii.gz"
Subepemask="${output_regionmask_dir}/Subepe_mask.nii.gz"
WMCSFmask="${output_regionmask_dir}/WMandCSF_mask.nii.gz"
GMCSFmask="${output_regionmask_dir}/GMandCSF_mask.nii.gz"
GMWMmask="${output_regionmask_dir}/GMandWM_mask.nii.gz"
NotGMmask="${output_regionmask_dir}/NotGM_mask.nii.gz"

echo "    Copying of Files & Calculating Masks is Done! Now Calculating Relevant TimeSeries!"
##############################################################################################################
fslmeants -i ${funcfile} -o "NotGM_${funcfilename}_timeseries.txt" -m ${NotGMmask}
fslmeants -i ${funcfile} -o "WM_${funcfilename}_timeseries.txt" -m ${WMmask}
fslmeants -i ${funcfile} -o "CSF_${funcfilename}_timeseries.txt" -m ${CSFmask}
fslmeants -i ${funcfile} -o "InnerCSF_${funcfilename}_timeseries.txt" -m ${InnerCSFmask}
fslmeants -i ${funcfile} -o "Suscept_${funcfilename}_timeseries.txt" -m ${Susceptmask}
fslmeants -i ${funcfile} -o "Outbrain_${funcfilename}_timeseries.txt" -m ${Outbrainmask}
fslmeants -i ${funcfile} -o "OutbrainOnly_${funcfilename}_timeseries.txt" -m ${OutbrainOnlymask}
fslmeants -i ${funcfile} -o "Edge_${funcfilename}_timeseries.txt" -m ${Edgemask}
fslmeants -i ${funcfile} -o "WMCSF_${funcfilename}_timeseries.txt" -m ${WMCSFmask}
fslmeants -i ${funcfile} -o "Subepe_${funcfilename}_timeseries.txt" -m ${Subepemask}

echo "    Relevant TimeSeries Are Computed! Now Running MELODIC!"
##################### MELODIC ICA DECOMPOSITION SECTION ######################################################
# let's finally get to melodic IC decomposition
cd ${output_dir}

# Need to calculate TR - use awk

TR="$(awk '{print $2}' <<< $(echo "$(fslhd ${funcfile})" | grep "pixdim4"))"

# MELODIC SECTION
# if melodic folder does not exist, then rerun it
if [ ! -d "${mel_fol}" ]; then
  mkdir "${mel_fol}"
  echo "      Running: melodic --in="${funcfile}" --outdir=${mel_fol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
  melodic --in="${funcfile}" --outdir=${mel_fol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
elif [ ! -f "${mel_fol}/melodic_IC.nii.gz" ]; then
  echo "      Melodic file is missing melodic_IC.nii.gz. Rerunning!"
  rm -rf "${mel_fol}"
  mkdir "${mel_fol}"
  echo "      Running: melodic --in="${funcfile}" --outdir=${mel_fol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}"
  melodic --in="${funcfile}" --outdir=${mel_fol} --mask="${funcmask}" --Ostats --nobet --mmthresh=0.5 --report --tr=${TR}
else
  echo "      Melodic already complete at ${mel_fol}"
fi

# Merge the probability maps and threshold maps which may prove useful later
probmapnames="$(ls ${mel_fol}/stats/probmap_* | sort -V)"
fslmerge -t ${mel_fol}/ICprobabilities.nii.gz ${probmapnames}

# Calculate a 99% ICprobabilities mask
fslmaths "${mel_fol}/ICprobabilities.nii.gz" -thr 0.99 -bin "${mel_fol}/ICprobabilities_99percent.nii.gz"

thresholdnames="$(ls ${mel_fol}/stats/thresh_zstat* | sort -V)"
fslmerge -t ${mel_fol}/ICthresh_zstat.nii.gz ${thresholdnames}

cd ${output_dir}

echo "    Melodic is Complete! Now Calculating Smoothness and Probabilities"
#########################################################################################
## Calculate smoothness retention before and after 6mm gaussian smoothing of IC
# Get explained variance
cd ${output_dir}
if [ -d "ROIcalcs" ]
then
  rm -rf ROIcalcs
fi
mkdir "ROIcalcs"
ROIcalcfol="${output_dir}/ROIcalcs"

# grab explained variance percentage in case of use in weighting the Data
awk '{print $1}' ${mel_fol}/melodic_ICstats > ${ROIcalcfol}/IC_exp_variance.txt

# Calculate the higher probabilities to use
fslmaths ${mel_fol}/ICprobabilities.nii.gz -thr 0.95 "${ROIcalcfol}/highprob_tmp_prob.nii.gz"
fslmaths "${ROIcalcfol}/highprob_tmp_prob.nii.gz" -bin "${ROIcalcfol}/highprob_tmp_mask.nii.gz"

# calculate a nonthresholded version of the z stats, absolute valued
fslmaths ${mel_fol}/ICthresh_zstat.nii.gz -abs "${ROIcalcfol}/fullvolICA_tmp_nothresh.nii.gz"

# calculate a smoothed 6mm gauss version of nonthresholded too
fslmaths ${mel_fol}/ICthresh_zstat.nii.gz -s 6 -abs "${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz"

# Calculate what you need form these to compare before and after smoothing to get a smoothing retention
fslstats -t "${ROIcalcfol}/fullvolICA_tmp_nothresh.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_nothresh_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_nothresh_ICnumvoxels.txt

fslstats -t "${ROIcalcfol}/fullvolICA_tmp_smoothed_nothresh.nii.gz" -M -V > ${ROIcalcfol}/curr_tmp_calc.txt
awk '{print $1}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_smoothed_nothresh_ICmean.txt
awk '{print $2}' ${ROIcalcfol}/curr_tmp_calc.txt > ${ROIcalcfol}/fullvolume_smoothed_nothresh_ICnumvoxels.txt

#########################################################################################
# Now set up brain network template for easy comparisons
echo "    Smoothing Checks are Complete! Now Setting Up Brain Network Template"
# check if the resampled template already exists, if not, then make it
template_dir="${cicada_home_dir}/templates"
if [ ! -d "${template_dir}" ]; then
  mkdir "${template_dir}"
fi

if [ ! -f "${template_dir}/network_template_${task_name}.nii.gz" ]; then
  echo "    Template Resampled Not Found, Resampling one now at ${template_dir}/network_template_${task_name}.nii.gz ."
  flirt -ref "${funcmask}" -in "${brain_network_template}" -out "${template_dir}/network_template_${task_name}.nii.gz" -usesqform -applyxfm -interp nearestneighbour
fi
# OK, now we have a resampled template that we can use, should have 1-7 values
# 1: Medial Visual, 2: Sensory Motor, 3: Dorsal Attention, 4: Ventral Attention, 5: FrontoParietal, 6: Default Mode Network, 7: Subcortical
network_template="${template_dir}/network_template_${task_name}.nii.gz"
#########################################################################################
# NOW, do all the relevant calculations!
echo "    Brain Network Template Calculated. Now finally do all the relevant calculations!"
calcfile="${ROIcalcfol}/highprob_tmp_prob.nii.gz"

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
rm -f ${output_dir}/*_tmp_*.nii.gz
rm -f ${ROIcalcfol}/*_tmp_*.nii.gz

echo "  Script 1 is Complete!"
echo
echo
###########################################################################################################
) 2>&1 | tee ${cicada_logfile}
