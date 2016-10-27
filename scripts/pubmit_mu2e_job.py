#!/bin/env python

import os
import subprocess
import argparse

class Submitter:

#------------------------------------------------------------------------------
    def __init__(self):
        self._parser = argparse.ArgumentParser();
        self._parser.add_argument("-c","--config"   , dest="fclFile"         , default="none" , help="FCL file name");
        self._parser.add_argument("-d","--dsid"     , dest="dsid"            , default="none" , help="input dataset ID");
        self._parser.add_argument("-f","--file"     , dest="fileName"        , default="none" , help="define filename");
        self._parser.add_argument("-F","--Files"    , dest="nFilesPerSegment", default=    -1 , help="N(files) per segment");
        self._parser.add_argument("-g","--grid"     , dest="gridCluster"     , default="local", help="GRID cluster");
        self._parser.add_argument("-n","--nevents"  , dest="nEvents"         , default=    -1 , help="N events per job");
        self._parser.add_argument("-N","--nseg"     , dest="nJobSegments"    , default=    -1 , help="N job segments");
        self._parser.add_argument("-o","--odir"     , dest="outputDir"       , default="none" , help="output directory");
        self._parser.add_argument("-s","--inputFile", dest="inputFile"       , default="none" , help="input file");
        self._parser.add_argument("-S","--inputDset", dest="inputDataset"    , default="none" , help="input dataset");
        self._parser.add_argument("-v","--verbose  ", dest="verbose"         , default="no"   , help="verbose mode");
        self._parser.add_argument("-x","--envVar"   , dest="envVar"          , default="no"   , help="env variables");

#        p = subprocess.Popen(['which','mu2e'],stdout=subprocess.PIPE,stderr=subprocess.PIPE);
#        out,err = p.communicate()
#        self._exefile = out.strip();

        self._exefile = subprocess.check_output('which mu2e',shell=True).strip();

#------------------------------------------------------------------------------
    def ParseCommandLine(self):
        self._args = self._parser.parse_args();
#------------------------------------------------------------------------------
# some parameters require additional parsing.
# 1. JOB is defined by the name of .fcl file
#------------------------------------------------------------------------------
        self._job = os.path.basename(self._args.fclFile).strip('.fcl')
        print 'JOB: %s'%self._job;

        self._dataset_dir = '/grid/fermiapp/mu2e/personal/murat/'+self._args.dsid;
# if the directory above exists, take a filelist from there

        self._dataset = self._dataset_dir;
        if (os.path.exists(dataset_dir)):
            self._dataset = dataset_dir+'/'+self._args.dsid+'.filelist';

#------------------------------------------------------------------------------
    def Execute(self):

        cmd = self._exefile

        print 'Execute: %s'%cmd


#------------------------------------------------------------------------------
# interactive submission
#------------------------------------------------------------------------------
if (__name__ == '__main__'):

    submitter = Submitter();
    submitter.ParseCommandLine();
    submitter.Execute();

#    subprocess.call(['./test.sh']);
