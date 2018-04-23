.PHONY : clean, default, publish

TEX_FILES   = limbdark
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out ent
LATEX       = pdflatex -interaction=nonstopmode -halt-on-error
SUFF        = pdf
RM_TMP      = $(foreach d, ${TEX_FILES}, rm -rf $(foreach suff, ${TMP_SUFFS}, ${d}.${suff}))
CHECK_RERUN = grep Rerun limbdark.log
FIGURES     = $(patsubst %.py,%.pdf,$(wildcard figures/*.py))
SOURCES     = $(wildcard ../src/*.[ch])

default: limbdark.pdf

limbdark.pdf: limbdark.tex stylez.tex $(FIGURES) $(SOURCES)
	${LATEX} limbdark.tex
	bibtex limbdark
	( ${CHECK_RERUN} && ${LATEX} limbdark.tex ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} limbdark.tex ) || echo "Done."

clean:
	$(RM_TMP)
	rm -rf limbdarkNotes.bib
	rm -rf $(FIGURES)

figures/%.pdf: figures/%.py
	cd $(<D);python $(<F)

publish:
	cd ..; ./publish.sh
