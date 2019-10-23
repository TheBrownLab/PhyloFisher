#!/bin/sh

javac $(pwd)/jap/*.java ;
javac $(pwd)/BMGE.java ;
if [ ! -e $(pwd)/BMGE.class ]; then exit; fi;

mkdir $(pwd)/JARBUILDER/
mkdir $(pwd)/JARBUILDER/jap/ ;
cp $(pwd)/jap/*.class $(pwd)/JARBUILDER/jap/ ;
cp -r $(pwd)/Jama/ $(pwd)/JARBUILDER/ ;
cp $(pwd)/*.class $(pwd)/JARBUILDER/ ;

mkdir $(pwd)/JARBUILDER/META-INF/ ;
echo "Main-Class: BMGE" > $(pwd)/JARBUILDER/META-INF/MANIFEST.MF ;
echo "Class-Path: ./Jama/*.class ./jap/*.class" >> $(pwd)/JARBUILDER/META-INF/MANIFEST.MF ;

cd ./JARBUILDER/ ;
jar cvfm BMGE.jar ./META-INF/MANIFEST.MF . ;
cd ../ ;

cp JARBUILDER/BMGE.jar . ;

rm -r $(pwd)/JARBUILDER/ ; 
rm $(pwd)/jap/*.class $(pwd)/*.class ;
