# Make errors fatal
set -e

# Get Travis CPU speed
echo "Travis CPU clock speed:"
lscpu | grep "MHz"

# Miniconda
export PATH="$HOME/miniconda-cache/bin:$PATH"

# Check if the conda command exists, and if not,
# download and install miniconda
if ! command -v conda > /dev/null; then

    # Conda
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b -p $HOME/miniconda-cache -u;
    conda config --add channels conda-forge;
    conda config --set always_yes yes;
    conda update --all;
    conda create --yes -n test python=$PYTHON_VERSION
    source activate test
    conda install tectonic;
    conda install -c conda-forge numpy=$NUMPY_VERSION scipy matplotlib setuptools pytest pytest-cov pip

    # Install exoplanet
    pip install astropy
    pip install exoplanet

    # Install batman
    pip install batman-package

    # Install a stable version of pytransit
    pushd $HOME
        git clone https://github.com/hpparvi/pytransit.git
        cd pytransit
        git checkout 7ce6eb238a64b29dc2001ff8b61311342820dfec # last stable commit for gimenez model
        python setup.py config_fc --fcompiler=gnu95 --opt="-Ofast" --f90flags="-cpp -fopenmp -march=native" build install
    popd

    # Install the dev version of starry (29-04-2019)
    pip install pillow
    pip install ipython
    pip install pybind11
    pushd $HOME
        git clone https://github.com/rodluger/starry.git
        cd starry
        git checkout c38c984593dd9e02f2dadc8ee11119852781b36e
        STARRY_BITSUM=515 python setup.py develop
    popd

fi

# Display some info
conda info -a

# Attempt to resolve issues with SSL certificate expiring for purl.org:
# https://tectonic.newton.cx/t/how-to-use-tectonic-if-you-can-t-access-purl-org/44
mkdir -p $HOME/.config/Tectonic
cat > $HOME/.config/Tectonic/config.toml << EOL
[[default_bundles]]
url = "http://purl.org/net/pkgwpub/tectonic-default"
EOL

# EXPERIMENTAL: Also installing texlive
# to resolve issues w/ Julia & LaTeX
texlive_year="2018"
sudo apt-get -qq update
export PATH=/tmp/texlive/bin/x86_64-linux:$PATH
if ! command -v pdflatex > /dev/null; then
    echo "Texlive not installed"
    echo "Downloading texlive and installing"
    wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
    tar -xzf install-tl-unx.tar.gz
    ./install-tl-*/install-tl --profile=./utilities/texlive.profile
    echo "Finished install TexLive"
fi
echo "Now updating TexLive"
tlmgr option -- autobackup 0
tlmgr update --self --all --no-auto-install
echo "Finished updating TexLive"
