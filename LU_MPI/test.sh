#!/bin/bash

if [ -n "$1" ]
then
	fl=1
	echo "------------------ONLY RESIDUAL AND TIME------------------"
else
	fl=0
	echo "-----------------------ALL RESULTS------------------------"
fi
	
dir="./SLU/small_tests"

data_s4="\
a.txt \
a20.txt \
b.txt \
"

data_s6="\
c.txt \
d.txt \
e.txt \
f.txt \
"



for (( proc = 1; proc <= 4; proc++ ))
do \
	for data in ${data_s4}
	do \
		data_file="${dir}/${data}"
		echo "====================== ${data_file} ======================"
		result=`mpirun -np ${proc} ./a.out 4 3 4 0 ${data_file} | tr -d '\0'`\
	
		if [ $fl -eq 1 ]
		then
			residual=`echo "${result}" | grep -a 'residual'`
		else
			echo "${result}"
		fi
		echo "------------------------------------------------------------"
	done

	for data in ${data_s6}
	do \
		data_file="${dir}/${data}"
		echo "====================== ${data_file} ======================"
		result=`mpirun -np ${proc} ./a.out 6 3 6 0 ${data_file} | tr -d '\0'`\
	
		if [ $fl -eq 1 ]
		then
			residual=`echo "${result}" | grep -a 'residual'`
			echo "${residual}"
		else
			echo "${result}"
		fi
		echo "------------------------------------------------------------"
	done
done
   






for (( proc = 1; proc <= 4; proc++ ))
do \
	for (( s = 1; s <= 4; s++ ))
	do \
		for (( n = 3; n <= 30; n++ ))
		do \
			for (( m = 3; m <= 30; m+=3 ))
			do \
					
				echo "================= n = ${n} m = ${m} s = ${s} ================= "
				result=`mpirun -np ${proc} ./a.out ${n} ${m} 10 ${s}| tr -d '\0'`\
					
				if [ $fl -eq 1 ]
				then	
					residual=`echo "${result}" | grep -a 'residual'`					
					echo "${residual}"
				else
					echo "${result}"
				fi					
			done
		done
	done
done






