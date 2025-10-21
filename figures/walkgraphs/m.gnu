# --- Terminal Setup ---
set title "Levyjevi pobegin in sprehodi pri različnih $\mu$"
if(ARG1 eq "png") {
    set terminal pngcairo size 1600,1200 enhanced font 'Verdana,10'
    set output 'm.png'
}
else {
    set terminal epslatex size 14cm,10cm
    set output 'm.tex'
}

# Use pm3d for color coding lines by Z-value (Column 3)
set palette model RGB defined ( 0 'blue', 1 'cyan', 3 'green', 5 'yellow', 7 'red' )
set pm3d map

# --- Multiplot Setup ---
set multiplot title "Levyjevi pobegin in sprehodi pri različnih $\mu$"

# --- 1. Calculate the Global Min/Max T for the color bar range ---
stats 'f1.dat' using 3 nooutput
stats 'f3.dat' using 3 nooutput
stats 'w1.dat' using 3 nooutput
stats 'w3.dat' using 3 nooutput
T_MIN = STATS_min
T_MAX = STATS_max
stats 'f3.dat' using 3 nooutput; if (STATS_min < T_MIN) {T_MIN = STATS_min}; if (STATS_max > T_MAX) {T_MAX = STATS_max}
stats 'w1.dat' using 3 nooutput; if (STATS_min < T_MIN) {T_MIN = STATS_min}; if (STATS_max > T_MAX) {T_MAX = STATS_max}
stats 'w3.dat' using 3 nooutput; if (STATS_min < T_MIN) {T_MIN = STATS_min}; if (STATS_max > T_MAX) {T_MAX = STATS_max}
set cbrange [T_MIN:T_MAX]


# --- 2. Configure and Plot the Color Bar (The Dummy Plot Method) ---

# Set a very small, temporary size/origin for the dummy plot
set size 0.001, 0.001
set origin -10.0, -10.0

# Configure the color box appearance and position
set colorbox user origin 0.92, 0.1 size 0.03, 0.8  # (x_start, y_start) (width, height)
set cblabel "Čas" rotate by 270 offset 3, 0

# The dummy plot command to FORCE the color bar to render
# It must plot something that uses the color palette (like a line from one of the files)
# The output is tiny and hidden, but it draws the color box.
plot 'f1.dat' using 1:2:3 with lines lc palette notitle

# --- 3. Plot the 4 Main Panels ---

# Restore the original plot size and origin settings
set size 1.0, 1.0
set origin 0.0, 0.0

# Turn OFF the color box for the individual panels (it's already placed)
unset colorbox

# --- Panel 1: f1.dat (Top Left) ---
set origin 0.0, 0.5   # x_start, y_start
set size 0.45, 0.45  # width, height (0.45 * 2 = 0.90 width)
set title "Pobeg, $\mu=1$"
set xlabel "x"
set ylabel "y"
plot 'f1.dat' using 1:2:3 with linespoints lc palette notitle


# --- Panel 2: f3.dat (Top Right) ---
set origin 0.45, 0.5
set size 0.45, 0.45
set title "Pobeg, $\mu=3$"
set xlabel "x"
unset ylabel # Remove Y label for cleaner look
plot 'f3.dat' using 1:2:3 with linespoints lc palette notitle


# --- Panel 3: w1.dat (Bottom Left) ---
set origin 0.0, 0.0
set size 0.45, 0.45
set title "Sprehod, $\mu=1$"
set xlabel "x"
set ylabel "y"
plot 'w1.dat' using 1:2:3 with linespoints lc palette notitle


# --- Panel 4: w3.dat (Bottom Right) ---
set origin 0.45, 0.0
set size 0.45, 0.45
set title "Sprehod, $\mu=3$"
set xlabel "x"
unset ylabel
plot 'w3.dat' using 1:2:3 with linespoints lc palette notitle


# --- Finalization ---

unset multiplot
pause -1
