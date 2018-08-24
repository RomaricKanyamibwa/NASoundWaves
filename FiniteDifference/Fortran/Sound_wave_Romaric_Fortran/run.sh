#!/bin/bash

if [ $# -ne 4 ]; then
    echo "usage : $0 [Nx] [Nt] [NumImages] [order]"
    exit 1
fi

Nx=$1
Nt=$2
NumImages=$3
Order=$4

echo "cleaning simu_impliciteV1/"
rm -rf simu_impliciteV1/*
# echo "Creating simu_impliciteV1/"
# mkdir -vp simu_implicite

Step=$((Nt/ NumImages))
echo "Step:$Step"

make
echo "./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages --order $Order"
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=false
./bin/SoundWaves -x $Nx -t $Nt --nimages $NumImages --order $Order

echo "gnuplot -c $load \"gnuplot_main.gp\" $Nx $Nt $Step $Order"

gnuplot -c $load "gnuplot_main.gp" $Nx $Nt $Step $Order


echo "*********************************"
echo "*      Leading Order system     *"
echo "*********************************"

echo "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"
xdg-open "Generated_files/Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"

if [ $Order -ge 1 ]; then
    echo "*********************************"
    echo "*      First Order system       *"
    echo "*********************************"
    echo "Generated_files/FOrder_Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"
    xdg-open "Generated_files/FOrder_Pressure_evolution_Nx_"$Nx"_Nt_"$Nt".gif"
fi
