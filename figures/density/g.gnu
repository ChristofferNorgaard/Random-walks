set terminal pngcairo size 800,600
set output 'density_heatmap.png'

# Create 2D bins
set dgrid3d 50,50 qnorm 2
set pm3d map
set view map
set palette rgbformulae 33,13,10

splot 'w15.dat' using 1:2:(1) with pm3d notitle
