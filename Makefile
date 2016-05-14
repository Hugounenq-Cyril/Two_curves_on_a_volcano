explicit-isogenies.pdf: explicit-isogenies.tex refs.bib \
benchmarks/graphe-101.pdf benchmarks/graphe-101-149-269.pdf
	latexmk -pdf explicit-isogenies.tex

benchmarks/101.pdf: benchmarks/101.dat benchmarks/101.gp
	cd benchmarks && gnuplot 101.gp

benchmarks/graphe-101.pdf: benchmarks/test-script-101-bis.tsv benchmarks/graphe-101.gp
	cd benchmarks && gnuplot graphe-101.gp

benchmarks/graphe-101-149-269.pdf: benchmarks/test-script-101-bis.tsv \
benchmarks/test-script-149-bis.tsv benchmarks/test-script-269-bis.tsv \
benchmarks/test-script-521-bis.tsv benchmarks/test-script-1033-bis.tsv \
benchmarks/test-script-62bits-bis.tsv benchmarks/test-script-252bits-bis.tsv \
benchmarks/graphe-101-149-269.gp
	cd benchmarks && gnuplot graphe-101-149-269.gp

