#!/bin/bash
set -eu -o pipefail
echo ${PREFIX}
mkdir -p ${PREFIX}/bin
make
cp dist_est ${PREFIX}/bin