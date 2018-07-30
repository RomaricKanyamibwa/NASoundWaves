# first animation 
# by Anida Khizar and Boussad MERHANE

###############################################################################
#     Generate animations using the generated files                           #
###############################################################################

reset
# set linetype 1 lc rgb "dark-violet" with linespoints
set term gif animate
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
# plot "simu_impliciteV1/out0_PH0.dat" using 1:2 title "Fluid pressure";
# plot "simu_impliciteV1/out1_PH0.dat" using 1:2 title "Fluid pressure";
while(i <= n){
	plot "simu_impliciteV1/out".i."_PH0.dat" using 1:2 title "Fluid pressure" with lines ,"simu_impliciteV1/out".i."_PH0.dat" using 1:3 title "Fluid velocity" with lines; 
	i = i+step;
}

set output
