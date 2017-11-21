#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

printenv

cd $VIRTUAL_ENV/tests

# Build
export PYTHONPATH=..
echo $PYTHONPATH

export PATH=../bin:$PATH
echo $PATH

# unit tests
python Test.py -v
