#!/opt/local/bin/gnuplot
#

# set terminal svg size 350,262 fname 'Verdana' fsize 10
# set output 'introduction.svg'

reset
# png
#set terminal pngcairo size 1024,768 enhanced font 'Verdana,9'
#set output 'log.png'

# svg
set terminal svg size 1024,768 fname 'Verdana, Helvetica, Arial, sans-serif' \
fsize '9' rounded dashed
set output 'log.svg'

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# color definitions
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green

set key bottom right

set xlabel 'log(t)'
set ylabel 'log(Psi(x, t))'

plot 'plot.dat' u 1:2 t 'Numerical solution' w lp ls 1, \
	""	u 1:3 t 'Accurate solution' w lp ls 2
