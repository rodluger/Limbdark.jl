echo "TESTING PYTRANSIT BEGIN"

pushd .ci
python -c "from test_pytransit import test; test(order=2, N=10, number=1)"
python -c "from test_pytransit import test; test(order=5, N=10, number=1)"
python -c "from test_pytransit import test; test(order=15, N=10, number=1)"

python -c "from test_pytransit import test; test(order=2, N=100, number=1)"
python -c "from test_pytransit import test; test(order=5, N=100, number=1)"
python -c "from test_pytransit import test; test(order=15, N=100, number=1)"

python -c "from test_pytransit import test; test(order=2, N=10000, number=1)"
python -c "from test_pytransit import test; test(order=5, N=10000, number=1)"
python -c "from test_pytransit import test; test(order=15, N=10000, number=1)"

python -c "from test_pytransit import test; test(order=2, N=10, number=10)"
python -c "from test_pytransit import test; test(order=5, N=10, number=10)"
python -c "from test_pytransit import test; test(order=15, N=10, number=10)"

python -c "from test_pytransit import test; test(order=2, N=10000, number=10)"
python -c "from test_pytransit import test; test(order=5, N=10000, number=10)"
python -c "from test_pytransit import test; test(order=15, N=10000, number=10)"

popd

echo "TESTING PYTRANSIT END"