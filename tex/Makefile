.PHONY : clean, default, publish

TEX_FILES   = limbdark
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out ent
LATEX       = pdflatex -interaction=nonstopmode -halt-on-error
SUFF        = pdf
RM_TMP      = $(foreach d, ${TEX_FILES}, rm -rf $(foreach suff, ${TMP_SUFFS}, ${d}.${suff}))
CHECK_RERUN = grep Rerun limbdark.log
PYFIGURES   = $(patsubst %.py,%.pdf,$(wildcard figures/python/*.py))
JLFIGURES   = $(patsubst %.jl,%.pdf,$(wildcard figures/julia/*.jl))

default: limbdark.pdf

limbdark.pdf: limbdark.tex stylez.tex $(JLFIGURES) $(PYFIGURES)
	${LATEX} limbdark.tex
	bibtex limbdark
	( ${CHECK_RERUN} && ${LATEX} limbdark.tex ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} limbdark.tex ) || echo "Done."

clean:
	$(RM_TMP)
	rm -rf limbdarkNotes.bib
	rm -rf $(PYFIGURES)
	rm -rf $(JLFIGURES)

figures/julia/%.pdf: figures/julia/%.jl
	cd $(<D);julia $(<F)

figures/python/%.pdf: figures/python/%.py
	cd $(<D);python $(<F)
