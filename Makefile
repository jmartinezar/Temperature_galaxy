temperature_galaxy.pdf: temperature_galaxy.tex 1.pdf 2.pdf
	pdflatex $<

1.pdf: plot.py include.py
	python3 $<

2.pdf: plot.py include.py
	python3 $<
