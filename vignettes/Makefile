ldDesign.pdf: ldDesign.Rnw ldDesign.bib
	$(R_HOME)/bin/R CMD Sweave ldDesign.Rnw
	latex ldDesign.tex
	bibtex ldDesign
	latex ldDesign.tex
	latex ldDesign.tex
	dvips ldDesign.dvi -o ldDesign.ps
	ps2pdf ldDesign.ps
	rm *.tex *.dvi *.ps *.aux *.log *.toc *.bbl *.blg

clean:
	rm *.tex *.dvi *.ps *.aux *.log *.toc *.bbl *.blg *.pdf
