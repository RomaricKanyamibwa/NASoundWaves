###############################################################################
#     Generate animations using the  files generated by the python script     #
###############################################################################

reset
set term gif animate
set title font 'Helvetica,14'
set title "Sound waves in rarefied gas around smooth solid bodies Nx=2000" 
set xlabel 'x'	
set ylabel 'y'

set output "./Sound_wave_Nx_2000_Nt_1000_to_9000_by_2000.gif"
#set n = `echo $n`
i = 1
n=1000
set xrange [0:100] #until xmax=D
set yrange [-4:10]

while(i < n){
	system "clear"
	print "### GIF creation ###"
	print "	",i*100/n,"%"
	plot "P_Nx_2000_Nt_2000_img_".(2*i).".txt" using 1:2 title "Pressure, Nt=2000" with lines linecolor rgb "black","P_Nx_2000_Nt_4000_img_".(i*4).".txt" using 1:2 title "Pressure, Nt=4000" with lines linecolor rgb "blue","P_Nx_2000_Nt_8000_img_".(i*8).".txt" using 1:2 title "Pressure, Nt=8000" with lines linecolor rgb "red","P_Nx_2000_Nt_16000_img_".(i*16).".txt" using 1:2 title "Pressure, Nt=16000" with lines linecolor rgb "green","P_Nx_2000_Nt_32000_img_".(i*32).".txt" using 1:2 title "Pressure, Nt=32000" with lines linecolor rgb "pink";
	i=i+1	
}

set output