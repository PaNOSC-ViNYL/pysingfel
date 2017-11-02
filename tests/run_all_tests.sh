#!/bin/sh
# USAGE: ./run_all_tests.sh 2>&1 | tee log
# This prints all messages and errors to the terminal and writes them to a file named log.

# Check if radiationDamageMPI is callable.
RADIATION_DAMAGE_EXECUTABLE=radiationDamageMPI
command -v ${RADIATION_DAMAGE_EXECUTABLE} >/dev/null 2>&1 || { echo >&2 "$RADIATION_DAMAGE_EXECUTABLE} was not found. Please make sure it's in your PATH."; exit 1; }

# Run the test suite.
python Test.py $@
