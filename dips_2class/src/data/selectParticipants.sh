#!/bin/bash

phenfile="phenotype/Phenotypic_V1_0b_preprocessed1.csv"
origin="selection/abide_subjects_func_ok_dimartino_tr_2_subjects_dx_group_without_missing.txt"
input="selection/abide_subjects_func_ok_dimartino_tr_2_subjects_dx_group_without_missing_id.txt"
output="selection/abide_subjects_func_ok_dimartino_tr_2_subjects_dx_group_without_missing_table.csv"

# extract subject IDs
# echo "matlab -nodesktop -nosplash -nojvm -r 'extractID('${origin}', '${input}'); exit;'"
matlab -nodesktop -r "extractID('${origin}', '${input}'); exit;"

# # extract
# > $output
# while read -r line;
# do
#     awk -v subjid=${line} \
#         '{ if ( $3 == subjid && $104 == "OK" ) print }' \
#         FPAT='([^,]*)|("[^"]*")' \
#         < $phenfile \
#         >> $output
# done < $input
