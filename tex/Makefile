.PHONY : clean, default, publish

TEX_FILES   = limbdark
TMP_SUFFS   = pdf aux bbl blg log dvi ps eps out ent
LATEX       = pdflatex -interaction=nonstopmode -halt-on-error
SUFF        = pdf
RM_TMP      = $(foreach d, ${TEX_FILES}, rm -rf $(foreach suff, ${TMP_SUFFS}, ${d}.${suff}))
CHECK_RERUN = grep Rerun limbdark.log
PYFIGURES   = $(patsubst %.py,%.pdf,$(wildcard figures/python/*.py))
JLFIGURES   = $(patsubst %.jl,%.pdf,$(wildcard figures/julia/*.jl))
TECTONIC    = $(shell command -v tectonic >/dev/null && echo true || echo false )

default: limbdark.pdf

limbdark.pdf: limbdark.tex stylez.tex $(JLFIGURES) $(PYFIGURES)
	# Generate links to current git commit
	python gitlinks.py
	if [ "${TECTONIC}" = "true" ]; then\
		tectonic limbdark.tex --print --keep-logs;\
	else\
		${LATEX} ${TEXFILE};\
		bibtex limbdark;\
		( ${CHECK_RERUN} && ${LATEX} ${TEXFILE} ) || echo "Done.";\
		( ${CHECK_RERUN} && ${LATEX} ${TEXFILE} ) || echo "Done.";\
		( ${CHECK_RERUN} && ${LATEX} ${TEXFILE} ) || echo "Done.";\
	fi

clean:
	$(RM_TMP)
	rm -rf limbdarkNotes.bib
	rm -rf $(PYFIGURES)
	rm -rf $(JLFIGURES)

figures/julia/s2_residuals_MA2002.pdf: figures/julia/s2_residuals_MA2002.jl
	cd figures/julia;\
	f2py -c occultquad.f -m occultquad;\
	julia s2_residuals_MA2002.jl

figures/julia/%.pdf: figures/julia/%.jl
	cd $(<D);julia $(<F)

figures/python/%.pdf: figures/python/%.py
	cd $(<D);python $(<F)
