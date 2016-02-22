set terminal post eps color enhanced font "Times,16" size 21cm,6cm
set out "creat.eps"

#set format y ""
set key right top
set rmargin 0
set lmargin 0
#set size square

set size 0.425,1
set origin 0.1,0
set ylabel "seconds" offset 5
set xlabel "r^2" offset -15,1.5
#set xtics 0, 2000
#set ytics 0, 5000
plot 'result.dat' index 0 using 1:2 with line lt 1 lw 3 title "p=101"
unset ylabel
unset xlabel
