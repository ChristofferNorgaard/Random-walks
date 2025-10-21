### 2D density color and contour plot from external data file
datafile = 'w35.dat'  # default filename

# Check if file exists and load data
t0 = time(0.0)
set table $Data
    plot datafile u 1:2 w table
unset table
print sprintf("Data loading: %.3f sec",(t1=time(0.0))-t0,t0=t1)

stats $Data u 0 nooutput   # get number of datalines
Rows = STATS_records
print sprintf("Loaded %d data points", Rows)

stats $Data u 1 nooutput
xmin = STATS_min
xmax_val = STATS_max
xmax = (abs(xmin) > abs(xmax_val)) ? abs(xmin) : abs(xmax_val)

stats $Data u 2 nooutput
ymin = STATS_min
ymax_val = STATS_max
ymax = (abs(ymin) > abs(ymax_val)) ? abs(ymin) : abs(ymax_val)

max_radius = sqrt(xmax**2 + ymax**2)

R_percent = 10.0            # Percentage of max radius (adjust as needed: 2-10%)
R = max_radius * R_percent / 100.0

print sprintf("Data range: x=[%.3f, %.3f], y=[%.3f, %.3f]", xmin, xmax_val, ymin, ymax_val)
print sprintf("Max radius: %.3f, Using R = %.3f (%.1f%%)", max_radius, R, R_percent)

# for each datapoint: how many other datapoints are within radius R
distance(x0,y0,x1,y1) = sqrt((x1-x0)**2 + (y1-y0)**2)
set print $Density
    do for [i=1:Rows] {
        x0 = real(word($Data[i],1))
        y0 = real(word($Data[i],2))
        c  = 0
        stats $Data u (distance(x0,y0,$1,$2)<=R ? c=c+1 : 0) nooutput
        d = c / (pi * R**2)
        print sprintf("%g %g %g", x0, y0, d)
    }
set print
print sprintf("Density check: %.3f sec",(t1=time(0.0))-t0,t0=t1)

set size ratio -1
set palette rgb 33,13,10
set cblabel "Število točk na enoto ploščine"

if(ARG1 eq "png") {
    set terminal pngcairo size 1600,1200 enhanced font 'Verdana,10'
    set output 'w35.png'
}
else {
    set terminal epslatex size 14cm,10cm
    set output 'w35.tex'
}

set title "Levyjevi sprehodi v navadnem difuzivnem režimu $\mu=3.5$"
set xrange[-1000:1000]
set yrange[-1000:1000]
plot $Density u 1:2:3 w p pt 7 ps 0.5 lc palette z notitle
