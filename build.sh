#!/usr/bin/env bash


# copy fisher main_scripts to package
mkdir -p $PREFIX/PhyloFisher
cp -r $RECIPE_DIR/phylofisher $PREFIX/PhyloFisher/phylofisher
cp -r $RECIPE_DIR/lib $PREFIX/PhyloFisher/lib


# Makes MAIN scripts executable and creates hard links in bin/
main_scripts=("fisher"\
              "config"\
              "fishing_net"\
              "forest"\
              "forge"\
              "informant"\
              "lumberjack"\
              "purge"\
              "build_dataset"\
              "trimming"\
              "missing_data" \
              )
for script in ${main_scripts[@]}; do
  chmod u+x $PREFIX/PhyloFisher/phylofisher/$script.py
  ln $PREFIX/PhyloFisher/phylofisher/$script.py $PREFIX/bin/$script.py
  done

# Makes UTILITY scripts executable and creates hard links in bin/
utility_scripts=("no_gap_stops"\
                 "add_aln_length"\
                 "sub"\
                 "pre_trimal"\
                 "AgentEllie"\
                 "AgentEllie2"\
                 "len_filter"\
                 "len_filter2"\
                 "fast_tax_removal"\
                 "heteroevolving_sites"\
                 )
for script in ${utility_scripts[@]}; do
  chmod u+x $PREFIX/PhyloFisher/phylofisher/utilities/$script.py
  ln $PREFIX/PhyloFisher/phylofisher/utilities/$script.py $PREFIX/bin/$script.py
  done


ln $PREFIX/PhyloFisher/lib/Divvier/divvier $PREFIX/bin/divvier
ln $PREFIX/PhyloFisher/lib/prequal/prequal $PREFIX/bin/prequal
ln $PREFIX/PhyloFisher/lib/BMGE-1.12/BMGE $PREFIX/bin/BMGE