#! /bin/bash
clear
make
clear
count=1
#10 tests
while [ $count -le $1 ]
do 
	echo "Now is test ${count}"
	nvprof ./main 1
	count=`expr $count + 1`
done
