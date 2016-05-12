#!/bin/bash

# defaults

andrei_dsid='cd3-cosmic-g4s4-general'

# 'cd3-cosmic-g4s4-general' or similar

if [ .$1 != "." ] ; then andrei_dsid=$1 ; fi

if   [ $andrei_dsid == 'cd3-cosmic-g4s4-general' ] ; then dsid=c4ga5760 ;
elif [ $andrei_dsid == 'cd3-cosmic-g4s4-target1' ] ; then dsid=c41a5760 ;
elif [ $andrei_dsid == 'cd3-cosmic-g4s4-target2' ] ; then dsid=c42a5760 ;
elif [ $andrei_dsid == 'cd3-cosmic-g4s4-target3' ] ; then dsid=c43a5760 ;
elif [ $andrei_dsid == 'cd3-cosmic-g4s4-target4' ] ; then dsid=c44a5760 ;
fi

d0=/pnfs/mu2e/persistent/datasets/phy-sim/sim/mu2e/${andrei_dsid}/v542_v566_v574_v576/art

for d1 in `ls $d0` ; do 
#    echo $d1
    d=$d0/$d1
    for d2 in `ls $d0/$d1` ; do
	dd=$d0/$d1/$d2
	ls $dd/*.art
    done
done >| datasets/$dsid/$dsid.files

