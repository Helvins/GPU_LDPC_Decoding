#! /bin/bash
clear
make
clear
count=1
#10 tests
while [ $count -le $1 ]
do 
	echo "Now is test ${count}"
	./main 0 20
	count=`expr $count + 1`
done
