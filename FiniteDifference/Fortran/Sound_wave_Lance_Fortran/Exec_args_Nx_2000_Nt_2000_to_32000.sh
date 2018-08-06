#!/bin/bash

rm file_txt/*.txt

clear

gfortran -O3 -ftree-vectorize -ffast-math implicit_scheme_second_order_for_GIF_variable_Nt.f90 -o implicit_scheme_second_order_for_GIF_variable_Nt

./implicit_scheme_second_order_for_GIF_variable_Nt 2000 2000
./implicit_scheme_second_order_for_GIF_variable_Nt 2000 4000
./implicit_scheme_second_order_for_GIF_variable_Nt 2000 8000
./implicit_scheme_second_order_for_GIF_variable_Nt 2000 16000
./implicit_scheme_second_order_for_GIF_variable_Nt 2000 32000

cd file_txt

gnuplot -c gnuplot_main_fixed_Nx_variable_Nt.gp
rm *.txt

xdg-open Sound_wave_Nx_fixed_Nt_variable_image.gif