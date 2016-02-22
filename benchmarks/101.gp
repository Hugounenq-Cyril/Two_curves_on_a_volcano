set terminal post eps color enhanced font "Times,16" size 8cm,4cm
set out "101.eps"

set key right bottom
set ylabel "seconds" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
plot '101.dat' using 1:($3+$4) with line lt 1 lw 3 title "diagonalization", \
     '101.dat' using 1:5 with line lt 2 lw 3 title "interpolation"
