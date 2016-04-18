set terminal post eps color enhanced font "Times,16" size 8cm,4cm
set out "graphe-101-521-1033.eps"

set key right bottom
set ylabel "seconds" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
plot 'test-script-101.tsv' using 1:4 with line lt 1 lw 3 title "p=101", \
     'test-script-521.tsv' using 1:4 with line lt 2 lw 3 title "p=521", \
     'test-script-1033.tsv' using 1:4 with line lt 3 lw 3 title "p=1033"
