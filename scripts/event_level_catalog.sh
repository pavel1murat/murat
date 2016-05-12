#!/bin/bash

dsid=cg4a5740
if [ .$1 != "." ] ; then dsid=$1 ; else dsid=cg4a5740 ; fi

catalog_dir=$PWD/datasets/$dsid

file_catalog=$catalog_dir/$dsid.files

if [ ! -d $catalog_dir/events ] ; then mkdir $catalog_dir/events ; fi

for file in `cat $file_catalog` ; do
# for file in /pnfs/mu2e/persistent/dsids/phy-sim/sim/mu2e/cd3-cosmic-g4s4-general/v542_v566_v574_v574/art/03/f1/sim.mu2e.cd3-cosmic-g4s4-general.v542_v566_v574_v574.001500_00400000.art ; do
#    echo $file
    of=`echo $file | awk -F / '{print $(NF-2)"_"$(NF-1)"_"$NF}'`

    mu2e -c murat/test/eventLister.fcl -s $file >| $catalog_dir/events/$of 2>&1 
done
