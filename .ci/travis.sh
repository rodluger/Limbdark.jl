# Miniconda (cached)
export PATH="$HOME/miniconda-cache/bin:$PATH"
if [ ! -f $HOME/miniconda-cache/bin/conda ]; then
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
      bash miniconda.sh -b -p $HOME/miniconda-cache -u;
      conda config --add channels conda-forge;
      conda config --set always_yes yes;
      conda update --all;
      conda create --yes -n test python=$PYTHON_VERSION
      conda activate test
      conda install tectonic;
      conda install -c conda-forge numpy=$NUMPY_VERSION scipy matplotlib setuptools pytest pytest-cov pip starry
fi

# Display some info
conda info -a

# Install batman
pip install batman-package

# Install pytransit
pushd $HOME
git clone https://github.com/hpparvi/pytransit.git
cd pytransit
python setup.py config_fc --fcompiler=gnu95 --opt="-Ofast" --f90flags="-cpp -fopenmp -march=native" build install
popd

# Install the dev version of starry
pip install pybind11
pushd $HOME
git clone https://github.com/rodluger/starry.git
cd starry
#git checkout 8cea59fa71acfcf4bd6265d6b21c38e82971436d
git checkout fb8462b762c03566168579a35d87abf6996df77c
STARRY_BITSUM=1 STARRY_KEEP_DFDU_AS_DFDG=1 python setup.py develop
popd

# Attempt to resolve issues with SSL certificate expiring for purl.org:
# https://tectonic.newton.cx/t/how-to-use-tectonic-if-you-can-t-access-purl-org/44
mkdir -p $HOME/.config/Tectonic
cat > $HOME/.config/Tectonic/config.toml << EOL
[[default_bundles]]
url = "http://purl.org/net/pkgwpub/tectonic-default"
EOL
