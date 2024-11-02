#!/bin/bash

for file in 1_YachidaS_2019 2_WirbelJ_2018 3_ZellerG_2014 4_VogtmannE_2016
do
	for a in yesedge noedge
	do
		if [ ${a} == yesedge ]
		then
			for pos in 5
			do
				for neg in 2
				do
					nohup python -u 0_CRC.py ${a} ${pos} ${neg} ${file} 50 &
				done
			done
		fi

		#if [ ${a} == noedge ]
		#then
		#	nohup python -u 0_CRC.py ${a} ${a} ${a} ${file} 50 &
		#fi
	done
done