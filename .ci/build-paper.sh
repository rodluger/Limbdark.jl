#!/bin/bash -x
# DEBUG set -e

# Only build the paper with Julia 0.7
if [ $TRAVIS_JULIA_VERSION == "0.7.0" ]
then

    # Generate the Julia figures
    # https://github.com/JuliaPy/PyPlot.jl/issues/317#issuecomment-337348563
    echo "Generating julia figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/julia
	f2py -c occultquad.f -m occultquad
    for f in *.jl; do
        echo "Running $f..."
        julia "$f"
    done

    # Generate the Python figures
    echo "Generating python figures..."
    cd $TRAVIS_BUILD_DIR/tex/figures/python
    for f in *.py; do
        echo "Running $f..."
        python "$f"
    done

    # Build the paper with tectonic
    cd $TRAVIS_BUILD_DIR/tex
    python gitlinks.py
	tectonic limbdark.tex

    # Force push the paper to GitHub
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
