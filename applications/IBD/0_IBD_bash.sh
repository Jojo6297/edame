#!/bin/bash

for file in 1_HMP_2019_ibdmdb 2_VilaAV_2018 3_HallAB_2017 4_NielsenHB_2014
do
	for a in yesedge noedge
	do
		if [ ${a} == yesedge ]
		then
			for pos in 4 5 6
			do
				for neg in 2
				do
					nohup python -u 0_IBD.py ${a} ${pos} ${neg} ${file} 50 > nohup_IBD_p${pos}_n${neg}_${file}.out &
				done
			done
		fi

		if [ ${a} == noedge ]
		then
			nohup python -u 0_IBD.py ${a} ${a} ${a} ${file} 50 > nohup_IBD_p${a}_n${a}_${file}.out &
		fi
	done
done
