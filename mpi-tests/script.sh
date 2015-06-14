#!/bin/bash

trap "exit" INT
for file in dense/*-in.gph 
do
    base=${file%-*}
    echo $base
    mpirun -np 4 ../Debug/mpi $file out 1>/dev/null
    sort out > outs
    touch ${base}.out
    sort ${base}.out > bases
    diff outs bases
done
rm out outs bases
