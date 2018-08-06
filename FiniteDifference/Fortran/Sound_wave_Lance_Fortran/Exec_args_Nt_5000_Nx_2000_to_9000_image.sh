#!/bin/bash

rm file_txt/*.txt

clear

gfortran -O3 -ftree-vectorize -ffast-math implicit_scheme_second_order_for_GIF_variable_Nx_image.f90 -o implicit_scheme_second_order_for_GIF_variable_Nx_image

./implicit_scheme_second_order_for_GIF_variable_Nx_image 2000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx_image 3000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx_image 5000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx_image 7000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx_image 9000 5000

cd file_txt

gnuplot -c gnuplot_main_fixed_Nt_variable_Nx_periode_image.gp

rm *.txt
