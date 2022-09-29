#!/bin/sh

# to copy over files from functional folders to CADICA folders

currsubjids=(102 103 105 106 107 108 109 112 114 115 116 117 121 122 125 126 128 129 130 132 134 135 136 138 139 140 142 143 144 145 146 147 149 153 156 157 158 160 161 164 165 168 169 171 172 173 174 176 178 179 181 184 185 186 187)
home="/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME"
preproc="${home}/Preproc_ICA_rest"
derivatives="${preproc}/derivatives"
fmriprepdata="${derivatives}/fmriprep"
taskid="rest"
CADICAfuncs="/Users/keithdodd/CADICA" # your folder containing the scripts and such
sessids=(01 02) # session numbering 01 02
TR="2" # repetition time of scan in seconds


for l in "${currsubjids[@]}"
do
  echo "Running for Subject ${l}"
  for j in "${sessids[@]}"
  do
  echo "Running for session ${j}"
  # go to current fmriprep subject folder
  currsubjfol="${fmriprepdata}/sub-${l}"
  cd ${currsubjfol}

  # find functional folder to copy from
  funcfol="${currsubjfol}/ses-${j}/func"

  # find CADICA folder to copy to
  CADICAfol="${derivatives}/CADICA"
  CADICAsubfol="${CADICAfol}/sub-${l}" # main subject folder to work within
  sessiondir="${CADICAsubfol}/ses-${j}"

  rm ${sessiondir}/confounds_timeseries
  # now copy over what you need/want
  cp "${funcfol}/sub-${l}_ses-${j}_task-${taskid}_desc-confounds_timeseries.tsv" ${sessiondir}/confounds_timeseries.csv

  done
done
