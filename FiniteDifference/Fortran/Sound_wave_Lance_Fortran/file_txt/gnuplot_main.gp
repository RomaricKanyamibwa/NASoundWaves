###############################################################################
#     Generate animations using the  files generated by the python script     #
###############################################################################

reset
set term gif animate
set title font 'Helvetica,14'
set title "Sound waves in rarefied gas around smooth solid bodies\nNx = ".ARG1." and Nt = ".ARG2 
set xlabel 'x'	
set ylabel 'y'

set output "./Sound_wave_Nx_".ARG1."_Nt_".ARG2.".gif"
#set n = `echo $n`
i = 1
n=ARG2
set xrange [0:100] #until xmax=D
set yrange [-4:8]

while(i < n){
	system "clear"
	print "### GIF creation ###"
	print "Nx = ",ARG1," and Nt = ",ARG2
	print "	",i*100/n,"%"
	if (i<10){
		plot "img_000".i.".txt" using 1:2 title "Pressure" with lines, "img_000".i.".txt" using 1:3 title "Velocity" with lines; 
	}
	else{
		if (i<100){
			plot "img_00".i.".txt" using 1:2 title "Pressure" with lines, "img_00".i.".txt" using 1:3 title "Velocity" with lines; 
		}
		else{
			if (i<1000){
				plot "img_0".i.".txt" using 1:2 title "Pressure" with lines, "img_0".i.".txt" using 1:3 title "Velocity" with lines; 
			}
			else{
				plot "img_".i.".txt" using 1:2 title "Pressure" with lines, "img_".i.".txt" using 1:3 title "Velocity" with lines; 
			}
		}
	}

	if (n<101){
		i = i + 1;
	}
	else{
		i = i+(ARG2/1000+1);
	}
	
}

set output
