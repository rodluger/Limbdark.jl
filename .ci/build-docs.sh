#!/bin/bash
#set -e

# Only make the docs with current Julia version
if [ $TRAVIS_JULIA_VERSION != "0.7.0" ]
then

    # Make the docs
    cd $TRAVIS_BUILD_DIR/docs
    julia make.jl

    # Begin
    branch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
    echo "Building docs from ${branch} branch..."

    # Get git hash
    rev=$(git rev-parse --short HEAD)

    # Go to the html build
    cd $TRAVIS_BUILD_DIR/docs/build/

    # Commit & force push back
    git init
    touch .nojekyll
    git add -f .nojekyll
    git add -f *
    git -c user.name='travis' -c user.email='travis' commit -m "rebuild gh-pages at ${rev}"
    git push -q -f https://$GITHUB_USER:$GITHUB_API_KEY@github.com/$TRAVIS_REPO_SLUG HEAD:gh-pages >/dev/null 2>&1

    # Return to the top level
    cd $TRAVIS_BUILD_DIR

fi