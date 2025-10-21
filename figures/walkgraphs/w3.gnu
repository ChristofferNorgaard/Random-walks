# --- Configuration ---
data_file = "w3.dat"

# Set the title and axis labels
set title "2D Path Visualization (Color-Coded by Time)"
set xlabel "X-coordinate (Column 1)"
set ylabel "Y-coordinate (Column 2)"
set grid

# Set aspect ratio to be square (often helpful for path visualizations)
set size square

# --- Color Coding (Time) ---

# Enable the color box to show the time scale
set colorbox
set cblabel "Time (t)" offset 2,0

# Define a color palette. 'pm3d' is implicitly used when 'lc palette' is called.
# Here we use the standard 'jet' like palette (often visually effective)
# set palette model RGB defined ( 0 'blue', 1 'cyan', 3 'green', 5 'yellow', 7 'red' )
set palette model RGB defined ( 0 'blue', 1 'cyan', 3 'green', 5 'yellow', 7 'red' )

# Automatically set the color range (cbrange) based on the min/max of column 3
set cbrange [*:*] 


# --- Plot Command ---

# Plot the data:
# 1. 'data_file' specifies the input file.
# 2. 'using 1:2:3' specifies columns: X (1), Y (2), and Color (3).
# 3. 'with linespoints' draws a line connecting the points, and markers on them.
# 4. 'lc palette' tells gnuplot to determine the color ('lc' = linecolor) from the Z-value (column 3) 
#    using the currently set palette.
plot data_file using 1:2:3 with linespoints lc palette title "Path Progression"

# Keep the plot window open (for interactive sessions)
pause -1
