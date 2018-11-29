# Miniconda (cached)
export PATH="$HOME/miniconda/bin:$PATH"
if ! command -v conda > /dev/null; then
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
      bash miniconda.sh -b -p $HOME/miniconda -u;
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

# Install some dependencies
pip install batman-package

pushd $HOME
git clone https://github.com/hpparvi/pytransit.git
cd pytransit
python setup.py config_fc --fcompiler=gnu95 --opt="-Ofast" --f90flags="-cpp -fopenmp -march=native" build
python setup.py install
popd

# Attempt to resolve issues with SSL certificate expiring for purl.org:
# https://tectonic.newton.cx/t/how-to-use-tectonic-if-you-can-t-access-purl-org/44
mkdir -p $HOME/.config/Tectonic
cat > $HOME/.config/Tectonic/config.toml << EOL
[[default_bundles]]
url = "http://purl.org/net/pkgwpub/tectonic-default"
EOL
