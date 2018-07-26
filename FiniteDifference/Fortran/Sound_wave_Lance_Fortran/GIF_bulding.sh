#!/bin/bash

for i in `seq 1000 1000 5000 `; do
	for j in `seq 1000 1000 10000 `; do
		./Exec_args_Nx_then_Nt.sh $i $j;
	done
done