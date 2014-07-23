default:
	pdflatex paper.tex
	
all: paper.pdf

TEX=$(shell find ./ -name "*.tex")

paper.pdf: *.sty *.cls *.bib $(TEX)
	pdflatex -interaction=nonstopmode paper.tex
	bibtex paper
	pdflatex -interaction=nonstopmode paper.tex
	pdflatex -interaction=nonstopmode paper.tex
#	rm -f *.out *.log *.cb *.cb2 *.blg *.bbl *.aux *.toc
	rm -f *.blg *.log *.out

.PHONY: clean default
clean: 
	rm -f paper.pdf

full: clean all

todo:
	@echo `grep -R TODO * | wc -l` todos left
	@grep -R TODO *
