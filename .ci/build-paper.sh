#!/bin/bash -x
set -e

# Are there changes in the tex directory?
if git diff --name-only $TRAVIS_COMMIT_RANGE | grep 'tex/'
then

    # Generate the Julia figures
    # https://github.com/JuliaPy/PyPlot.jl/issues/317#issuecomment-337348563
    echo "Generating julia figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/julia
    for f in *.jl; do
        echo "Running $f..."
        LD_PRELOAD=${HOME}/.julia/v0.6/Conda/deps/usr/lib/libz.so julia "$f"
    done

    # Generate the Python figures
    echo "Generating python figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/python
    for f in *.py; do
        echo "Running $f..."
        python "$f"
    done

    # Install texlive
    sudo apt-get -qq update && sudo apt-get install -y --no-install-recommends texlive-full

    # Build the paper
    cd $TRAVIS_BUILD_DIR/tex
	pdflatex -interaction=nonstopmode -halt-on-error limbdark.tex
	bibtex limbdark
	pdflatex -interaction=nonstopmode -halt-on-error limbdark.tex
	pdflatex -interaction=nonstopmode -halt-on-error limbdark.tex
    pdflatex -interaction=nonstopmode -halt-on-error limbdark.tex

    # Force push the paper to GitHub
    # NOTE: This is hacky. I couldn't get the old aproach to work for some reason.
    cd $HOME
    mkdir tmp && cd tmp
    git init
    git checkout --orphan $TRAVIS_BRANCH-pdf
    mkdir tex
    cp $TRAVIS_BUILD_DIR/tex/limbdark.pdf tex/
    git add -f tex/limbdark.pdf
    git -c user.name='travis' -c user.email='travis' commit -m "building the paper"
    git push -q -f https://$GITHUB_USER:$GITHUB_API_KEY@github.com/$TRAVIS_REPO_SLUG $TRAVIS_BRANCH-pdf

fi
