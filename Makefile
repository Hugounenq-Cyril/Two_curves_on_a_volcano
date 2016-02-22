explicit-isogenies.pdf: explicit-isogenies.tex refs.bib benchmarks/101.eps
	latexmk -pdf explicit-isogenies.tex

benchmarks/101.eps: benchmarks/101.dat benchmarks/101.gp
	cd benchmarks && gnuplot 101.gp
