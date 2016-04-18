set terminal post eps color enhanced font "Times,16" size 8cm,4cm
set out "graphe-101-149-269.eps"

set key right bottom
set ylabel "seconds" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
plot 'test-script-101.tsv' using 1:4 with line lt 1 lw 3 title "p=101", \
     'test-script-149.tsv' using 1:4 with line lt 2 lw 3 title "p=149", \
     'test-script-269.tsv' using 1:4 with line lt 3 lw 3 title "p=269"
