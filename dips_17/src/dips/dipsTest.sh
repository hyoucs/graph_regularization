#!/bin/bash

j=0
nb_concurrent_processes=10
for idxFold in {1..10}
do

	echo "    ====== START: DIPS cv fold: ${idxFold} ======"
	matlab -nodesktop -r "dipsTest(${idxFold}); exit;" &
	((++j == nb_concurrent_processes)) && { j=0; wait; }
	echo "    ====== END: DIPS cv fold: ${idxFold} ======"

done