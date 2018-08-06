#!/bin/bash

rm file_txt/*.txt

clear
gfortran -O3 -ftree-vectorize -ffast-math implicit_scheme_second_order_for_GIF_variable_Nt_image.f90 -o implicit_scheme_second_order_for_GIF_variable_Nt_image

echo "Loading Nx = 2000 / Nt = 16000"
./implicit_scheme_second_order_for_GIF_variable_Nt_image 2000 16000
echo "Done Nx = 2000 / Nt = 16000"
clear
echo "Loading Nx = 2000 / Nt = 24000"
./implicit_scheme_second_order_for_GIF_variable_Nt_image 2000 24000
echo "Done Nx = 2000 / Nt = 24000"
clear
echo "Loading Nx = 2000 / Nt = 32000"
./implicit_scheme_second_order_for_GIF_variable_Nt_image 2000 32000
echo "Done Nx = 2000 / Nt = 32000"
clear
echo "Loading Nx = 2000 / Nt = 40000"
./implicit_scheme_second_order_for_GIF_variable_Nt_image 2000 40000
echo "Done Nx = 2000 / Nt = 40000"
clear
echo "Loading Nx = 2000 / Nt = 48000"
./implicit_scheme_second_order_for_GIF_variable_Nt_image 2000 48000
echo "Done Nx = 2000 / Nt = 48000"
clear

cd file_txt

gnuplot -c gnuplot_main_fixed_Nx_variable_Nt_periode_image.gp

rm *.txt
