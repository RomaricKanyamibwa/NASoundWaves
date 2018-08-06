#!/bin/bash

clear

gfortran -O3 -ftree-vectorize -ffast-math temperature_implicit_scheme.f90 -o temperature_implicit_scheme

./temperature_implicit_scheme $1 $2

cd file_txt

gnuplot -c temperature_gnuplot_main.gp $1 $2

#rm *.txt

xdg-open Temperature_sound_wave_Nx_$1_Nt_$2.gif