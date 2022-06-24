MAIN = using_drtmle_obstudies

all : tex pdf clean

tex:
	Rscript -e "rmarkdown::render('$(MAIN).Rmd')"
	sed -i '1s/^/\\include{preamble}\n /' $(MAIN).tex

pdf :
	pdflatex $(MAIN).tex
	bibtex $(MAIN)
	pdflatex $(MAIN).tex
	pdflatex $(MAIN).tex

clean :
	rm -f $(addprefix $(MAIN),.blg .fff .log .mw .out .aux .bbl)
	rm -f $(addprefix preamble,.aux)
