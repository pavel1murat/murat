#!/usr/bin/env python
#------------------------------------------------------------------------------
#  use --caf=local to submit a local job instead of going to CAF
#  so far comment out SAM Dbserver part....
#------------------------------------------------------------------------------
import string, getopt, os, sys, time, glob, commands, subprocess

#------------------------------------------------------------------------------
class TestSubmit:

    def __init__(self):
        self.fProject          = ''
        self.tarfile           = ''
        self.fNJobs            = 1
        self.fCaf               = 'local'
        self.fCalibPass        = ''
        self.fScript           = None 
        self.fSetupScript      = None
        self.fVerbose          = 0
	self.fDoit             = 1;
        self.lockfile          = None
        self.fOutputDir        = None;
        self.procType          = 'long'
        self.fUser             = 'cdfprd'
        self.email             = ''
        self.fJobNumber        = 1    ;   # for debugging
        self.fInputFile        = '';
        self.fWebTarball       = ''
  
#------------------------------------------------------------------------------
    def Time(self,Format='%Y/%m/%d %H:%M:%S'):
        return time.strftime(Format, time.localtime(time.time()))

    def Print(self,Name,Message):
        pid = str(os.getpid())
        message = self.Time()+'[cS::'+Name+'('+pid+')] '+Message
        print message

    def PrintCmdError(self,name,cmd,output,rc):
        self.Print(name,'ERROR: Command %s failed with rc=%s, output=%s' % (cmd,rc,output))        

    def ParseParameters(self):
        _name = 'ParseParameters'

        self.Print(_name,'arguments, sys.argv:[1:1] ');
        print sys.argv


        try:
            optlist, args = getopt.getopt(sys.argv[1:], '',
                                          ['project=', 'conf=', 'tarfile=',
                                           'samglia=', 'metalia=', 'size=', 'caf=',
                                           'station=', 'pass=', 'submit=', 'job=',
                                           'setup=', 'file=', 'doit=','print=',
                                           'verbose=','projsource=']) 

        except getopt.GetoptError, e:
            self.Print(_name,"%s %s" % ('argument error: ',sys.argv))
            sys.exit(1)


        for key, val in optlist:
            print  'args: ', key, val
            
            if key == '--project':
                self.fProject = val
            if key == '--job':
                self.fJobNumber = int(val)
            if key == '--file':
                self.fInputFile = val
            if key == '--verbose':
                self.fVerbose = int(val)
            if key == '--doit':
                self.fDoit = int(val)
            elif key == '--conf':
                self.DSConf = val
            elif key == '--tarfile':
                self.tarfile = val
            elif key == '--samglia':
                self.fSamgliaDir = val
            elif key == '--metalia':
                self.fMetaliaDir = val
            elif key == '--pass':
                self.fCalibPass = val
            elif key == '--size':
                self.fNJobs = int(val)
            elif key == '--caf':
                self.fCaf = val
            elif key == '--station':
                self.fSamStation = val
            elif key == '--submit':
                self.fScript = val
            elif key == '--setup':
                self.fSetupScript = val
            elif key == '--print':
                task = val
                if (task == 'projects'):
                    self.fVerbose = 1;
                    rc            = self.CheckProjects(self.fSamStation,0)
                    sys.exit(rc)
            elif key == '--projsource':
                self.fProjSource = val


        self.Print(_name, 'Done');

        return 0;

#------------------------------------------------------------------------------
    def Init(self):
        _name = 'Init'
        rc    = 0;


#------------------------------------------------------------------------------
# sequence of steps during the job submission is as follows:
#   1) start SAM project
#   2) submit the CAF job
#   3) declare job metadata using file of a dataset 'farm0i' as a container
#
# so if job metadata exists the SAM project should be already running....
#------------------------------------------------------------------------------
    def Run(self):

        _name = 'Run'
        rc    = 0;

#------------------------------------------------------------------------------
    def Finish(self):
        _name = 'Finish'
        
        self.Print(_name,'deleting lockfile %s before exit' % self.lockfile)
        os.remove(self.lockfile);

        return 0;



#------------------------------------------------------------------------------
def test_001(Project = 'CalibExe_1_bphysr'):
    x = TestSubmit();
    rc          = x.Init();

#    subprocess.
    

#------------------------------------------------------------------------------
#  Starting point, decide if need to proceed at all
#------------------------------------------------------------------------------
if __name__ == '__main__':
 
    print '\n\n  test_submit: %s ' % sys.argv
    
    rc = test_001();

    sys.exit(rc)
#------------------------------------------------------------------------------

