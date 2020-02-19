#!/usr/bin/env bash

mkdir $PREFIX/Prequal
cp -r $SRC_DIR/* $PREFIX/Prequal

ln $PREFIX/Prequal/prequal $PREFIX/bin/prequal
