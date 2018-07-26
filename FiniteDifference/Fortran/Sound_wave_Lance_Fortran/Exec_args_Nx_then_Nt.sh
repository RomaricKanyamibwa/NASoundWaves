#!/bin/bash

rm file_txt/*.txt

clear

gfortran -O3 -ftree-vectorize -ffast-math implicit_scheme_second_order.f90 -o implicit_scheme_second_order

./implicit_scheme_second_order $1 $2

cd file_txt

gnuplot -c gnuplot_main.gp $1 $2

rm *.txt

#xdg-open Sound_wave_Nx_$1_Nt_$2.gif