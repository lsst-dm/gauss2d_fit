#!/bin/sh
# This script ensures that pytest will succeed even if this executable is
# called on MacOS with SIP (which strips out DYLD_LIBRARY_PATH)
if [ "$LSST_LIBRARY_PATH" != "" ]; then
    export DYLD_LIBRARY_PATH=${LSST_LIBRARY_PATH}
fi
pytest tests "$@"
