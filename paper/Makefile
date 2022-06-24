MAIN = using_drtmle_obstudies

all : pdf clean

pdf:
	Rscript -e "rmarkdown::render('$(MAIN).Rmd')"

clean :
	rm -f $(addprefix $(MAIN),.blg .fff .log .mw .out .aux .bbl)
