#!/bin/bash

if [ $# -ne 3 ]; then
    echo "usage : $0 [Nx] [Nt] [NumImages]"
    exit 1
fi

Nx=$1
Nt=$2
NumImages=$3

echo "Deleting simu_impliciteV1/"
rm -rf simu_impliciteV1/*
mkdir simu_impliciteV1

Step=$((Nt/ NumImages))
echo "Step:$Step"

make
echo "./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages"

./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages

echo "gnuplot -c $load \"gnuplot_main.gp\" $Nx $Nt $Step"

gnuplot -c $load "gnuplot_main.gp" $Nx $Nt $Step

echo "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"
xdg-open "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"

