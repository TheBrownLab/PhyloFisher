#!/bin/bash
set -eu -o pipefail
mkdir -p ${PREFIX}/bin
make CC=$CC CPP=$CXX F77=$F77

cp mammal ${PREFIX}/bin
cp mammal-sigma ${PREFIX}/bin
cp mult-data ${PREFIX}/bin
cp mult-mix-lwt ${PREFIX}/bin
cp charfreq ${PREFIX}/bin
cp dgpe ${PREFIX}/bin