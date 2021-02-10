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

# Directories to process
pydirs="litebirdms"

# Run

for pyd in ${pydirs}; do
    find "${base}/${pyd}" -name "*.py" -exec ${blkexe} ${blkrun} '{}' \;
done
