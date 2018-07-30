#!/bin/bash

if [ $# -ne 3 ]; then
    echo "usage : $0 [Nx] [Nt] [NumImages]"
    exit 1
fi

Nx=$1
Nt=$2
NumImages=$3

Step=$((Nt/ NumImages))
echo "Step:$Step"

make
echo "./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages"

./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages

echo "gnuplot -c $load \"gnuplot_main.gp\" $Nx $Nt $Step"

gnuplot -c $load "gnuplot_main.gp" $Nx $Nt $Step

