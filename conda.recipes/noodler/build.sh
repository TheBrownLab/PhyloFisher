#!/usr/bin/env bash

mkdir ${PREFIX}/Noodler
mkdir ${PREFIX}/bin
cp -r ${SRC_DIR}/* ${PREFIX}/Noodler/

cat > ${PREFIX}/Noodler/Noodler.sh <<- EOM
#!/bin/bash

parent_path=\$( cd "\$(dirname "\${BASH_SOURCE[0]}")" ; pwd -P )
cd "\$parent_path"
index_html=../Noodler/index.html

firefox "\$index_html"

EOM

cp ${PREFIX}/Noodler/Noodler.sh ${PREFIX}/bin/Noodler
