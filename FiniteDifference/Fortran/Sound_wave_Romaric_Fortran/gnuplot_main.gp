# first animation 
# by Anida Khizar and Boussad MERHANE

###############################################################################
#     Generate animations using the generated files                           #
###############################################################################

reset
# set linetype 1 lc rgb "dark-violet" with linespoints
set term gif animate
set boxwidth 0.3 absolute
set style fill   solid 1.00 border lt -1
set grid nopolar
set grid xtics nomxtics ytics nomytics noztics nomztics nortics nomrtics \
 nox2tics nomx2tics noy2tics nomy2tics nocbtics nomcbtics
set grid layerdefault   lt 0 linecolor 0 linewidth 0.500,  lt 0 linecolor 0 linewidth 0.500
set key fixed right top vertical Right noreverse enhanced autotitle columnhead box lt black linewidth 1.000 dashtype solid
set key opaque
set style textbox opaque margins  1.0,  1.0 fc  bgnd border  lt -1 linewidth  1.0
set pointsize 2
set xtics border in scale 0,0 mirror norotate  autojustify
set xtics  norangelimit 
set xtics   ()
set ytics border in scale 0,0 mirror norotate  autojustify
set ztics border in scale 0,0 nomirror norotate  autojustify
set cbtics border in scale 0,0 mirror norotate  autojustify
set rtics axis in scale 0,0 nomirror norotate  autojustify
set title font 'Helvetica,14'
set title "Sound waves in Rarefied Gas: Fluid Pressure/velocity over time \nNx = ".ARG1." and Nt = ".ARG2 
set xlabel 'X'	
set ylabel 'Y'

set output "Generated_files/Pressure_evolution_Nx_".ARG1."_Nt_".ARG2.".gif"
#set n = `echo $n`
i = 0
n=ARG2
step=ARG3

print "Nx = ",ARG1,",Nt = ",ARG2 ," and Step=",ARG3
print "Generating gif image . . . . . ."

set xrange [-2:100] #until xmax=D
set yrange [-4:16]
# plot "simu_impliciteV1/out0_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:2 title "Fluid pressure";
# plot "simu_impliciteV1/out1_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:2 title "Fluid pressure";
while(i <= n){
#out0_Nx5000_Nt5000_PH0.dat
	plot "simu_impliciteV1/out".i."_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:2 title "Fluid pressure" with lines ,"simu_impliciteV1/out".i."_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:3 title "Fluid velocity" with lines , "simu_impliciteV1/out".i."_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:4 title "Fluid WH0" with lines, "simu_impliciteV1/out".i."_Nx".ARG1."_Nt".ARG2."_PH0.dat" using 1:5 title "Temperature" with lines; 
	i = i+step;
}

set output
