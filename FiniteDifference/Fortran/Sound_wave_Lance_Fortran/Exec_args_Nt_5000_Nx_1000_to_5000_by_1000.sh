#!/bin/bash

rm file_txt/*.txt

clear

gfortran -O3 -ftree-vectorize -ffast-math implicit_scheme_second_order_for_GIF_variable_Nx.f90 -o implicit_scheme_second_order_for_GIF_variable_Nx

./implicit_scheme_second_order_for_GIF_variable_Nx 1000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx 2000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx 3000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx 4000 5000
./implicit_scheme_second_order_for_GIF_variable_Nx 5000 5000

cd file_txt

gnuplot -c gnuplot_main_fixed_Nt_variable_Nx.gp

rm *.txt

xdg-open Sound_wave_Nt_5000_Nx_1000_to_5000_by_1000.gif