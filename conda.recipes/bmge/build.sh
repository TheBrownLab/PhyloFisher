#!/usr/bin/env bash

mkdir $PREFIX/BMGE-1.12
cp -r $SRC_DIR/* $PREFIX/

ln $PREFIX/BMGE $PREFIX/bin/BMGE
