#!/bin/bash
#------------------------------------------------------------------------------
# book-keeping for TS3 scan: 
# run on 5 files, N(POT) = 1e6 (200,000 POT, ot 20 subruns, per file)
# input: subruns 0000-0019 from concatenated Stage1 dataset: cd3-beam-cs1-mubeam.0506a
#------------------------------------------------------------------------------
export list=""
export doit=""

init_list() {
    x=$1

    if [[ ".$x" == "." || ".$x" == ".." ]] ; then 
	x='  70_40_0100 70_40_0200 70_40_0300 70_40_0400 70_40_0500 70_40_1000 '
	x=$x'80_40_0100 80_40_0200 80_40_0300 80_40_0400 80_40_0500 80_40_1000 '
	x=$x'80_50_0100 80_50_0200 80_50_0300 80_50_0400 80_50_0500 80_50_1000 '
	x=$x'90_50_0100 90_50_0200 90_50_0300 90_50_0400 90_50_0500 90_50_1000 '
	x=$x'90_60_0100 90_60_0200 90_60_0300 90_60_0400 90_60_0500 90_60_1000 '
    fi
#------------------------------------------------------------------------------
# x='80_40_0400 80_40_0500 80_40_1000 '
# x='70_40_0200' # debug
#------------------------------------------------------------------------------
    export list=$x
}

#------------------------------------------------------------------------------
make_fcl() {

    init_list $1
    doit=$2

    echo list=$list

    for x in $list ; do
	template_fcl=murat/test/ts3_tooth/g4s23/beam_g4s23_$x.fcl
	echo $template_fcl
	
	generate_fcl --desc=beam_g4s23_$x --dsconf=v0 \
            --inputs=/mu2e/app/users/murat/v6_1_4_prof/datasets/cd3-beam-cs1-mubeam.0506a/filesets/cd3-beam-cs1-mubeam.0506a.001002_00000000_00000099 \
	    --merge=1 $template_fcl

	mv seeds.murat.beam_g4s23* 000
	od=/mu2e/app/users/murat/fcl/ts3_tooth/beam_g4s23/$x
	mv 000 $od
	ls $od/*.fcl >| $od/aaa_fcl_list.txt
    done
}

#------------------------------------------------------------------------------
submit_jobs() {

    init_list $1
    doit=$2

    echo list=$list

    for x in $list ; do
	cmd="mu2eprodsys --setup=$PWD/setup.sh --wfpro=beam_g4s23_$x --fcllist=/mu2e/app/users/murat/fcl/ts3_tooth/beam_g4s23/$x/aaa_fcl_list.txt --dsconf=v622 --expected-lifetime=8h --memory=3GB"
	echo $cmd
	if [ ".$doit" != "." ] ; then $cmd ; fi
    done
}

#-----------------------------------------------------------------------------
handle_output() {

    init_list $1
    doit=$2

    if [ ."$doit" == "." ] ; then prefix=echo ; else prefix='' ; fi

    echo list=$list

    for x in $list ; do
	grid_dir="/pnfs/mu2e/scratch/users/murat/workflow/beam_g4s23_$x/outstage"
	odir="/mu2e/data/users/murat/datasets/ts3_tooth/v3/g4s23/$x"

	if [ ! -d $odir ] ; then mkdir -p $odir ; fi
#------------------------------------------------------------------------------
# always copy log files to keep the record
#------------------------------------------------------------------------------
	mkdir $odir/log
	cp $grid_dir/*/*/*/*.log  $odir/log/.
#------------------------------------------------------------------------------
# make file catalogs for different streams
#------------------------------------------------------------------------------
	mkdir $odir/catalog
	for stream in "mothers" "ootstops" "tgtstops" "truncated" ; do
	    ls $grid_dir/*/*/*/*.art | grep $stream >| $odir/catalog/g4s23_${stream}_$x
	done
    done    
}

#------------------------------------------------------------------------------
# grep on the path name, not module name
#------------------------------------------------------------------------------
calculate_muon_stops() {

    init_list $1
    doit=$2

    for x in $list ; do 
	odir="/mu2e/data/users/murat/datasets/ts3_tooth/v3/g4s23/$x"
	n=`cat $odir/log/*.log | grep TrigReport | grep tgtFilter | \
	    awk 'BEGIN{n=0} {if ($NF == "tgtFilter") n+=$5} END{print n}'`
	echo $x "  " $n
    done
}