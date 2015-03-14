#!/usr/bin/env python
#------------------------------------------------------------------------------
#  so far comment out SAM Dbserver part....
#------------------------------------------------------------------------------
import string, getopt, os, sys, time, glob, commands

#------------------------------------------------------------------------------
class Mu2eSubmit:

    def __init__(self):
        self.tarfile           = ''
        self.fCaf              = 'local'
        self.fJob              = None
        self.fDsid             = None
        self.fScript           = None 
        self.fSetupScript      = None
        self.fVerbose          = 0
	self.fDoit             = 1;
        self.fOutputDir        = None;
        self.fUser             = 'murat'
        self.fSegmentNumber    = 1   ;   # for debugging
        self.fInputFile        = '';
  
#------------------------------------------------------------------------------
    def Time(self,Format='%Y/%m/%d %H:%M:%S'):
        return time.strftime(Format, time.localtime(time.time()))


#------------------------------------------------------------------------------
    def Print(self,Name,Message):
        print self.Time()+' [ Mu2eSubmit::'+Name+' ] '+Message

#------------------------------------------------------------------------------
    def ParseParameters(self):
        _name = 'ParseParameters'

        self.Print(_name,'arguments, sys.argv:[1:1] ');
        print sys.argv


        try:
            optlist, args = getopt.getopt(sys.argv[1:], '',
                                          ['dsid=', 'segment=', 'job=',
                                           'file=', 'doit='   , 'print=',
                                           'verbose='])

        except getopt.GetoptError, e:
            self.Print(_name,"%s %s" % ('argument error: ',sys.argv))
            sys.exit(1)


        for key, val in optlist:
            print  'args: ', key, val
            
            if key == '--dsid':
                self.fDsid = val
            if key == '--job':
                self.fJob = val;
            if key == '--segment':
                self.fSegmentNumber = int(val)
            if key == '--file':
                self.fInputFile = val
            if key == '--verbose':
                self.fVerbose = int(val)
            if key == '--doit':
                self.fDoit = int(val)


        datetime.date dt;
        
        self.fOutputDir = "./results/"+dt.today().strftime("%Y-%m-%d");
        self.Print(_name, 'Done');

        return 0;

#------------------------------------------------------------------------------
    def Init(self):
        _name = 'Init'
        rc    = 0;

	if (self.fVerbose != 0): self.Print(_name,'Entered')
 
        rc = self.ParseParameters();

        if (rc != 0) : return rc;

        return rc;

#------------------------------------------------------------------------------
    def GetDatasets(self):
        _name = 'GetDatasets'
        

        return 0;

#------------------------------------------------------------------------------
    def Run(self):

        _name = 'Run';
        
        rc    = 0;

        self.Print(_name,'Entered')

        logfile = self.fDataset+".log"
       
        self.GetDatasets();
        sub  = 'source /grid/fermiapp/products/mu2e/setupmu2e-art.sh; source ./setup.sh ; ';
        sub += "cd $OUTPUT_DIR ; if [ -f "+logfile+" ] ; then rm "+logfile
        sub += 'mu2e -c '+self.fFclFile;
        sub += ' -S '+self.fDataset;
        
        if (self.fNEvents):
            sub += ' -n '+self.fNEvents;

        if (self.fCaf == 'local') :
                
            # running locally
                
            cmd = sub+(' --job=%i' % self.fJobNumber)
        else:
            x = 0;
#------------------------------------------------------------------------------
#  actually submit the job
#------------------------------------------------------------------------------
        cmd += " > $logfile 2>> $logfile "

        print "cmd = \n", cmd;

        self.Print(_name,'self.fDoit = %i, cmd = %s' % (self.fDoit,cmd))
        if (self.fDoit != 0):
            #------------------------------------------------------------------------------
            rc, output = commands.getstatusoutput(cmd)
            self.Print(_name, 'output= %s' % output)

                      
        return int(rc);


#------------------------------------------------------------------------------
    def Finish(self):
        _name = 'Finish'
        
        self.Print(_name,'deleting lockfile %s before exit' % self.lockfile)
        os.remove(self.lockfile);

        return 0;

    
#------------------------------------------------------------------------------
def test1(Project = 'CalibExe_1_bphysr'):
    js = Mu2eSubmit();
    js.fProject = Project;
    js.fVerbose = 1;
    js.fDoit    = 0;
    rc          = js.Init();
    
#------------------------------------------------------------------------------
#  Starting point, decide if need to proceed at all
#------------------------------------------------------------------------------
if __name__ == '__main__':
 
    print 'starting Mu2eSubmit: %s ' % sys.argv

    
    js = Mu2eSubmit();
    
    rc = js.Init();
    #--------------------------------------------------------------------------
    # run only if rc = 0,
    #--------------------------------------------------------------------------
    if ((rc == 0) and (js.fResourcesMissing == 0)):
        rc = js.Run();

    if (rc >= 0) :
        rc = js.Finish();

    sys.exit(rc)
#------------------------------------------------------------------------------

