#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

printenv

export TEST_DIR=$TRAVIS_BUILD_DIR/tests
cd $TEST_DIR

# Build
export PYTHONPATH=$TRAVIS_BUILD_DIR:$PYTHONPATH
echo $PYTHONPATH

export PATH=$TRAVIS_BUILD_DIR/bin:$PATH
echo $PATH

# unit tests
python Test.py -v
