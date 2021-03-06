###############################################################################
#     Generate animations using the  files generated by the python script     #
###############################################################################

reset
set term png
set title font 'Helvetica,14'
set xlabel 'x'	
set ylabel 'y'
#set n = `echo $n`
i = 0
n=10

while(i < n){

	set title "Sound waves in rarefied gas around smooth solid bodies\nperiod: pi/2+".i."*pi/2"
	set xrange [0:5+i] #until xmax=D
	set yrange [-2:2]
	set output "./Sound_wave_Nt_fixed_Nx_variable_image_k=".i.".png"
	plot "P_Nx_2000_Nt_5000_img_".(5000*(1+i)/400).".txt" using 1:2 title "Pressure, Nx=2000" with lines linecolor rgb "black","P_Nx_3000_Nt_5000_img_".(5000*(1+i)/400).".txt" using 1:2 title "Pressure, Nx=3000" with lines linecolor rgb "blue","P_Nx_5000_Nt_5000_img_".(5000*(1+i)/400).".txt" using 1:2 title "Pressure, Nx=5000" with lines linecolor rgb "red","P_Nx_7000_Nt_5000_img_".(5000*(1+i)/400).".txt" using 1:2 title "Pressure, Nx=7000" with lines linecolor rgb "green","P_Nx_9000_Nt_5000_img_".(5000*(1+i)/400).".txt" using 1:2 title "Pressure, Nx=9000" with lines linecolor rgb "brown"
	
	set output

	i=i+1
}
