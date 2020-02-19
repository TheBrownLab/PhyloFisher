#!/usr/bin/env bash

mkdir $PREFIX/BMGE-1.12
cp -r $SRC_DIR/* $PREFIX/BMGE-1.12

ln $PREFIX/BMGE-1.12/BMGE $PREFIX/bin/BMGE
