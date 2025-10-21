# Plot for Ai function errors
# set terminal pngcairo size 1600,1200 enhanced font 'Verdana,10'
# set terminal epslatex size 14cm,10cm
if(ARG1 eq "png") {
    set terminal pngcairo size 1600,1200 enhanced font 'Verdana,10'
    set output 'max_asymp.png'
}
else {
    set terminal epslatex size 14cm,10cm
    set output 'max_asymp.tex'
}
set datafile separator ","

set title 'Airyjeva funkcija Ai: Primerjava napake' font 'Verdana,14'
set xlabel 'x'
set ylabel 'Absolutna napaka' 
set logscale y
set grid

set key outside below

plot 'err.dat' using 1:2 with lines linewidth 2 title 'Maclaurinova vrsta', \
     'err.dat' using 1:3 with lines linewidth 2 title 'Asimptotska vrsta'

# Plot for Bi function errors
set output 'ex_bi_error.tex'

set title 'Airyjeva funkcija Bi: Primerjava napake' font 'Verdana,14'

plot 'err.dat' using 1:4 with lines linewidth 2 title 'Maclaurinova vrsta', \
     'err.dat' using 1:5 with lines linewidth 2 title 'Asimptotska vrsta'
