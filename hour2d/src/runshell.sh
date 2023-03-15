#!/bin/bash

N=3
# 1 2 4 8 16 32 64 128 192 256 512 1024
block_sizes=(8 16 32 64 128 192 256)
Best_block=0
Best_runtime=1000
arg="$*"

echo -e "\n`date`">>result_log.txt
echo -e "-----------------------------------------------------------------">>result_log.txt
for loop in ${block_sizes[@]}
do
	tmp=0
	sum=0
	exec="$arg OP_BLOCK_SIZE=$loop"
	declare -a runtime_list
	
	echo $exec>>result_log.txt
	while(($tmp < $N))
	do
		echo "$tmp $exec"
		echo "`$exec`">result_temp.txt
		runtime=`sed -n '$,$p' result_temp.txt`
		runtime_list[$tmp]=$runtime
		sum=`echo "scale=3;$sum + $runtime" | bc`
		let "tmp++"
	done
	echo " ${runtime_list[@]}">>result_log.txt
    sum=`echo "scale=3;$sum / $N" | bc`
	echo "ave:$sum">>result_log.txt
	if [ `echo "$sum<$Best_runtime"|bc` -eq 1 ] ; then
		Best_runtime=$sum
		Best_block=$loop
	fi
done
echo "Best_runtime:${Best_runtime}     Best_block:${Best_block}"
echo "Best_runtime:${Best_runtime}     Best_block:${Best_block}">>result_log.txt
echo "end">>result_log.txt

