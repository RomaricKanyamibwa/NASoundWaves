#!/bin/bash

if [ $# -ne 5 ]; then
    echo "usage : $0 [Nx] [Nt1] [Nt2] [Nt3] [NumImages]"
    exit 1
fi

Nx=$1
Nt1=$2
Nt2=$3
Nt3=$4
NumImages=$5

echo "cleaning simu_impliciteV1/"
rm -rf simu_impliciteV1/*
echo "Creating simu_impliciteV1/"
mkdir -vp simu_implicite

# Step=$((Nt/ NumImages))
# echo "Step:$Step"

make
echo "Generating .dat files"

./bin/SoundWaves -x $Nx -t $Nt1 --nimages $NumImages
./bin/SoundWaves -x $Nx -t $Nt2 --nimages $NumImages
./bin/SoundWaves -x $Nx -t $Nt3 --nimages $NumImages

echo "gnuplot -c $load \"gnuplot_errorNt.gp\" $Nx $Nt1 $Nt2 $Nt3 $NumImages"

gnuplot -c $load "gnuplot_errorNt.gp" $Nx $Nt1 $Nt2 $Nt3 $NumImages

echo "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt_"$Nt1"-$Nt2-$Nt3.gif"
xdg-open "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt"$Nt1"-$Nt2-$Nt3.gif"

