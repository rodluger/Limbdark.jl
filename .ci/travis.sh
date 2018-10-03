# Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

# Conda Python
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda create --yes -n test python=$PYTHON_VERSION
conda activate test
conda install -c conda-forge numpy=$NUMPY_VERSION scipy matplotlib setuptools pytest pytest-cov pip starry

# Install starry [DISABLED: Let's use conda-forge]
# pip install git+https://github.com/rodluger/starry

# Install texlive
sudo apt-get -qq update && sudo apt-get install -y --no-install-recommends texlive-full
