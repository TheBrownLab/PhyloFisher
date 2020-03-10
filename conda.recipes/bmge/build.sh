#!/usr/bin/env bash
mkdir $PREFIX/BMGE-1.12
cp -r * $PREFIX/BMGE-1.12/
mkdir $PREFIX/bin
ln $PREFIX/BMGE-1.12/BMGE $PREFIX/bin/BMGE