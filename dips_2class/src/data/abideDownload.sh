#!/bin/bash

subjects="./selection/abide_subjects_func_ok_dimartino_tr_2_subjects_dx_group_without_missing_table.csv"
# subjects="selection/abide_subset_ok_subjects.csv"

pipeline=cpac
strategy=filt_noglobal
derivative=rois_cc200
ext=1D

mkdir -p download
while read -r subjline;
do
	subj_id=$(cut -f2 -d, <<< $subjline);
    file_id=$(awk '{ print $7 }' FPAT='([^,]*)|("[^"]*")' <<< $subjline)
    url=https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/${pipeline}/${strategy}/${derivative}/${file_id}_${derivative}.${ext}

    wget \
        -O download/${pipeline}_${strategy}_${derivative}_${file_id}_${derivative}.${ext} \
        $url
done < $subjects
