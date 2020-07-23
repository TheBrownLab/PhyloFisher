#!/usr/bin/env bash

mkdir ${PREFIX}/parasorter
mkdir ${PREFIX}/bin

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        cp -r ${SRC_DIR}/parasorter-linux-x64/* ${PREFIX}/parasorter/
elif [[ "$OSTYPE" == "darwin"* ]]; then
        cp -r ${SRC_DIR}/parasorter-darwin-x64/* ${PREFIX}/parasorter
fi

ln -s ${PREFIX}/parasorter/parasorter ${PREFIX}/bin/parasorter
