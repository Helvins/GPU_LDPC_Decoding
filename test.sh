#! /bin/bash
clear
make
clear
count=1
#10 tests
while [ $count -le $1 ]
do 
	echo "Now is test ${count}"
	nvprof --unified-memory-profiling off --metrics achieved_occupancy ./main 10 1
	count=`expr $count + 1`
done
