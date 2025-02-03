#!/bin/bash

for file in 1_LifeLinesDeep_2016 2_AsnicarF_2021 3_MehtaRS_2018 4_ZeeviD_2015
do
	for a in yesedge noedge
	do
		if [ ${a} == yesedge ]
		then
			for pos in 4 5 6
			do
				for neg in 2
				do
					nohup python -u 0_healthy.py ${a} ${pos} ${neg} ${file} 50 > nohup_healthy_p${pos}_n${neg}_${file}.out &
				done
			done
		fi

		if [ ${a} == noedge ]
		then
			nohup python -u 0_healthy.py ${a} ${a} ${a} ${file} 50 > nohup_healthy_p${a}_n${a}_${file}.out &
		fi
	done
done
