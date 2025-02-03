#!/bin/bash

for file in 1_MetaCardis_2020_a 2_QinJ_2012 3_LiJ_2014 4_KarlssonFH_2013
do
	for a in yesedge noedge
	do
		if [ ${a} == yesedge ]
		then
			for pos in 6
			do
				for neg in 2
				do
					nohup python -u 0_T2D.py ${a} ${pos} ${neg} ${file} 50 > nohup_T2D_p${pos}_n${neg}_${file}.out &
				done
			done
		fi

		if [ ${a} == noedge ]
		then
			nohup python -u 0_T2D.py ${a} ${a} ${a} ${file} 50 > nohup_T2D_p${a}_n${a}_${file}.out &
		fi
	done
done
