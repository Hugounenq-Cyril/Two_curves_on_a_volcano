set terminal post eps color enhanced font "Times,16" size 8cm,4cm
set out "graphe-101-149-269.eps"

set key right bottom
set ylabel "seconds" offset 5
set xlabel "r"
set logscale y 2
set logscale x 2
plot 'test-script-101-bis.tsv' using 1:4 with line lt 1 lw 3 title "p=101", \
     'test-script-149-bis.tsv' using 1:4 with line lt 2 lw 3 title "p=149", \
     'test-script-269-bis.tsv' using 1:4 with line lt 3 lw 3 title "p=269", \
     'test-script-521-bis.tsv' using 1:4 with line lt 4 lw 3 title "p=521", \
     'test-script-1033-bis.tsv' using 1:4 with line lt 5 lw 3 title "p=1033", \
     'test-script-62bits-bis.tsv' using 1:4 with line lt 6 lw 3 title "p=62bits", \
     'test-script-252bits-bis.tsv' using 1:4 with line lt 7 lw 3 title "p=252bits", 
