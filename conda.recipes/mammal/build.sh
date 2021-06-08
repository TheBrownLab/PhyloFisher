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
cp *dat ${PREFIX}/bin

mkdir -p ${PREFIX}/etc/conda/activate.d

cat > ${PREFIX}/etc/conda/activate.d/${PKG_NAME}_activate.sh <<EOF
#!/bin/bash

if [ -d "~/.mammal" ]; then rm -Rf ~/.mammal;
  mkdir ~/.mammal
fi

cp ${PREFIX}/bin/*.dat ~/.mammal
EOF
