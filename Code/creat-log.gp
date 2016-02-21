set terminal post eps color enhanced font "Times,16" size 21cm,6cm
set out "creat.eps"

#set logscale y 2
set key right bottom
set rmargin 0
set lmargin 0
#set size square

set size 0.425,1
set origin 0.1,0
set ylabel "seconds" #offset 5
set xlabel "$r^2$" #offset -15,1.5
set xtics 0, 1
set ytics 0, 1
plot [0:10][0:10] 'result-log.dat' index 0 using 1:2 with line lt 1 lw 3 title "p=101",\
'result-log.dat' index 1 using 1:2 with line lt 2 lw 3 title "p=269",\
'result-log.dat' index 2 using 1:2 with line lt 3 lw 3 title "p=512",\
'result-log.dat' index 3 using 1:2 with line lt 4 lw 3 title "p=1033"
unset ylabel
unset xlabel
