echo "TESTING PYTRANSIT BEGIN"

pushd .ci
python -c "from test_pytransit import test; test(1, 2)"
python -c "from test_pytransit import test; test(100, 2)"
python -c "from test_pytransit import test; test(10000, 2)"

python -c "from test_pytransit import test; test(1, 5)"
python -c "from test_pytransit import test; test(100, 5)"
python -c "from test_pytransit import test; test(10000, 5)"

python -c "from test_pytransit import test; test(1, 10)"
python -c "from test_pytransit import test; test(100, 10)"
python -c "from test_pytransit import test; test(10000, 10)"

python -c "from test_pytransit import test; test(1, 15)"
python -c "from test_pytransit import test; test(100, 15)"
python -c "from test_pytransit import test; test(10000, 15)"
popd

echo "TESTING PYTRANSIT END"