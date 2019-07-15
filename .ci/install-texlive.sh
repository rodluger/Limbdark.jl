# Make errors fatal
set -e

# EXPERIMENTAL: Also installing texlive
# to resolve issues w/ Julia & LaTeX
texlive_year="2018"
sudo apt-get -qq update
export PATH=/tmp/texlive/bin/x86_64-linux:$PATH
if ! command -v pdflatex > /dev/null; then
    echo "Texlive not installed"
    echo "Downloading texlive and installing"
    wget --quiet http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
    tar -xzf install-tl-unx.tar.gz
    ./install-tl-*/install-tl --profile=./utilities/texlive.profile
    echo "Finished install TexLive"
fi
echo "Now updating TexLive"
tlmgr option -- autobackup 0
tlmgr update --self --all --no-auto-install
echo "Finished updating TexLive"
