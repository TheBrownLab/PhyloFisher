#!/usr/bin/env bash

mkdir $PREFIX/Divvier
cp -r $SRC_DIR/* $PREFIX/Divvier

ln $PREFIX/Divvier/divvier $PREFIX/bin/divvier
