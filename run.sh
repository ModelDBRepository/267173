#!/bin/bash

make 

mkdir data


function la_run {
	 echo ./lamodel $LAPARAMS   
	 ./lamodel $LAPARAMS   
}


./lamodel -V 

for run in {0..10}  ; do
	(
	LAPARAMS="  -P 2 -S 19$run  -s control_${run} -G  "
	la_run

	#LAPARAMS="  -P 2 -S 19$run  -s blocked_${run} -G -o nDA2Pyr=1000  -n  "
	LAPARAMS="  -P 2 -S 19$run  -s blocked_${run} -G -n "
	la_run
	) &

done

wait

python3 graphs.py

