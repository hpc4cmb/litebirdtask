#!/bin/bash
#
# Apply Python formatting with the "black" command line tool.
#

# Get the directory containing this script
pushd $(dirname $0) > /dev/null
base=$(pwd -P)
popd > /dev/null

# Check executables

blkexe=$(which black)
if [ "x${blkexe}" = "x" ]; then
    echo "Cannot find the \"black\" executable.  Is it in your PATH?"
    exit 1
fi

# Black runtime options
blkrun="-l 88"

# Black test options
blktest="--check"

# Note that the "+" argument to "find ... -exec" below passes all found files to the
# exec command in one go.  This works because black accept multiple files as arguments.

find "${base}/litebirdtask" -name "*.py" -and -not \
    -path '*/_version.py' -exec ${blkexe} ${blkrun} '{}' + &

# Wait for multiple commands to finish
wait
