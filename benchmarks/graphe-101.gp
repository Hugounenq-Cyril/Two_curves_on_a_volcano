set terminal pdf color enhanced font "Times,16" size 10cm,7cm
set out "graphe-101.pdf"

set key left top Left
set ylabel "seconds" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
set xtics 4
plot [4:4096] 'test-script-101-bis.tsv' using 1:($2+$3) with line lt 2 dt 2 lw 3 title "Algorithms 1+2", \
     'test-script-101-bis.tsv' using 1:4 with line lt 1 dt 1 lw 3 title "Interpolation"
