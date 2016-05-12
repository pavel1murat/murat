#!/bin/bash 


# file=/mu2e/app/users/murat/dev2_prof/datasets/cg4a5740/events/01_8d_sim.mu2e.cd3-cosmic-g4s4-general.v542_v566_v574_v574.001500_00244000.art

#------------------------------------------------------------------------------
help () {
    echo "call format: prepare_file_event_tables.sh DSID "
    echo "DSID is not defined, EXIT"
}

dsid=""
if [ .$1 != "." ] ; then dsid=$1 ; else dsid=cg4a5740 ; fi

if [ .$dsid == "." ] ; then 
    help 
    exit -1
fi

echo DSID=$dsid

catalog_dir=$PWD/datasets/$dsid

for file in `ls $catalog_dir/events/*.art` ; do
    fn=`cat $file | grep "Closed file" | awk '{print $NF}'`
#    echo $fn
    cat $file | grep "Event:" | awk -v f=$fn '{printf "%5i %7i %8i %s\n", $4,$6,$8,f}' 
done >| $catalog_dir/$dsid.events

