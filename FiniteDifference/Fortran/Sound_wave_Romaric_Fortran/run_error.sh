#!/bin/bash

if [ $# -ne 5 ]; then
    echo "usage : $0 [Nx] [Nx1] [Nx2] [Nt] [NumImages]"
    exit 1
fi

Nx=$1
Nx1=$2
Nx2=$3
Nt=$4
NumImages=$5

echo "cleaning simu_impliciteV1/"
rm -rf simu_impliciteV1/*
echo "Creating simu_impliciteV1/"
mkdir -vp simu_implicite

Step=$((Nt/ NumImages))
echo "Step:$Step"

make
echo "Generating .dat files"

./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages
./bin/SoundWaves -x $Nx1 -t $Nt --nimages $NumImages
./bin/SoundWaves -x $Nx2 -t $Nt --nimages $NumImages

echo "gnuplot -c $load \"gnuplot_error.gp\" $Nx $Nx1 $Nx2 $Nt $Step"

gnuplot -c $load "gnuplot_error.gp" $Nx $Nx1 $Nx2 $Nt $Step

echo "Generated_files/Pressure_evolution_Nx_"$Nx"-$Nx1-"$Nx2"_Nt_"$Nt".gif"
xdg-open "Generated_files/Pressure_evolution_Nx_"$Nx"-$Nx1-"$Nx2"_Nt_"$Nt".gif"

