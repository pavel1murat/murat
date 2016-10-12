#!/usr/bin/env python
#
#
# $Id: SConscript,v 1.7 2013/09/09 03:54:15 murat Exp $
# $Author: murat $
# $Date: 2013/09/09 03:54:15 $
#
import os, re, string, subprocess
Import('env')
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
x = subprocess.call('scripts/build_config',shell=True)
#------------------------------------------------------------------------------
# dictionary build rules
#------------------------------------------------------------------------------
def murat_gen_rootcint(source, target, env, for_signature):
#    print "\n>>> murat_gen_rootcint called:"
#    print ">>> source = ",len(source),str(source[0]), str(source[1])
#    print ">>> target = ",len(target),str(target[0])
#    print ">>> for_signature = ",for_signature

    class_include = str(source[1]);
    linkdef       = str(source[0]);

#    print "[stntuple_gen_rootcint] class_include = %s"%class_include
#    print "[stntuple_gen_rootcint] linkdef       = %s"%linkdef
    
    includes =   "-Iinclude -I"+os.environ['ART_INC'     ];
    includes = includes + " -I"+os.environ['CETLIB_INC'  ];
    includes = includes + " -I"+os.environ['CANVAS_INC'  ];
    includes = includes + " -I"+os.environ['FHICLCPP_INC'];

    dict        = str(target[0]);
    tmp_lib_dir = os.path.dirname(dict);                # env['MURAT_TMP_LIB_DIR']
    
    pcm_file = dict.replace(".cxx","_rdict.pcm");

#    print "dict:"+dict + "   pcm_file:"+pcm_file;
    lib_dir = os.environ['MU2E_BASE_RELEASE']+"/lib";

    cmd = 'if [ ! -d '+tmp_lib_dir+' ] ; then mkdir -p '+tmp_lib_dir+'; fi ;';
    cmd = cmd+"rootcint -f "+dict+" -c -DMU2E "+includes+" "+class_include+" "+linkdef+"; ";
    cmd = cmd+'if [ ! -d '+lib_dir+' ] ; then mkdir '+lib_dir+' ;';
    cmd = cmd+"mv "+pcm_file+" "+lib_dir+'/.'; 
#    print ">>> cmd = %s"%cmd
    return cmd

murat_my_rootcint = Builder(generator     = murat_gen_rootcint,
                            single_source = 0,
                            suffix        = '.o',
                            src_suffix    = '.h')

env.Append(BUILDERS = {'MuratRootCint' : murat_my_rootcint})


class murat_helper:
    """mu2e_helper: class to produce library names"""
#   This appears to behave like c++ static member and is initialized at class defintion time.
    sourceroot =  os.path.abspath('.')

    def __init__(self,debug = False):
        self._list_of_object_files = [];

        self.dd      = re.search('[^/]*/[^/]*$',env.Dir('.').abspath).group(0)
        self.dirname = os.path.dirname(self.dd);   # THIS
        self.subdir  = os.path.basename(self.dd);
        self.libname = self.dirname+'_'+self.subdir
        self.d1      = self.libname+'-shared';
        self.suffix  = ".hh"
        self._debug  = debug

        env['MURAT_TMP_LIB_DIR'] = 'tmp/src/'+self.d1

#
#   Accesor
#
    def base(self):
        return self.sourceroot

    def handle_dictionaries(self,suffix = ".hh"):
        self.suffix = suffix
#------------------------------------------------------------------------------
# generate dictionaries
#------------------------------------------------------------------------------
        list_of_linkdef_files = Glob(self.subdir+'/dict/*_linkdef.h', strings=True)
#        print "[Stntuple/obj] list_of_linkdef_files = \n",list_of_linkdef_files
        list_of_dict_files    = []

        for f in list_of_linkdef_files:
            linkdef       = string.split(str(f),'/');
            clname        = string.replace(linkdef[len(linkdef)-1],"_linkdef.h","");
            include       = self.subdir+'/'+clname+suffix;
            
            dict          = '#/tmp/src/'+self.d1+'/'+clname+'_dict.cxx';
            list_of_dict_files.append(dict);

            #    print "linkdef = ",linkdef
            #    print "include = ",include

            env.MuratRootCint(dict,[f,include])
#------------------------------------------------------------------------------
# compile dictionaries
#------------------------------------------------------------------------------
        list   = [];

        for dict in list_of_dict_files:
            obj_cxx_file = string.replace(dict,".cxx",".o");
            list.append(obj_cxx_file);
            self._list_of_object_files.append(obj_cxx_file);
        
            include  = string.replace(dict,".cxx",".h");
            env.SharedObject(obj_cxx_file,dict)

    def compile_fortran(self,list_of_f_files, skip_list = []):
        if (self._debug):
            print  self.dirname+"[build_libs]: list_of_f_files:"+self.subdir,list_of_f_files

        for f in list_of_f_files:
            if (not f in skip_list):
                o = '#tmp/src/'+self.d1+'/'+string.split(f,'.')[0]+'.o'
                self._list_of_object_files.append(o);
                env.SharedObject(o,f)

    def build_libs(self,list_of_cc_files, skip_list = [],libs = []):
        if (self._debug):
            print  self.dirname+"[build_libs]: list_of_cc_files:"+self.subdir,list_of_cc_files

        for cc in list_of_cc_files:
            if (not cc in skip_list):
                #    print ".cc file: "+cc
                o = '#/tmp/src/'+self.d1+'/'+string.split(cc,'.')[0]+'.o'
                self._list_of_object_files.append(o);
                env.SharedObject(o,cc)

        #        libs     = [ 'Stntuple_base' ]

        lib_name = os.environ['MU2E_BASE_RELEASE']+'/lib/'+self.libname+'.so';

        # print "Stntuple/obj list_of_obj_files:",list_of_obj_files

        env.SharedLibrary(lib_name,self._list_of_object_files,LIBS = [libs])
        # print ">>>>> EXIT Stntuple/obj"

    def build_modules(self,list_of_module_files, skip_list, libs = []):
        if (self._debug):
            print self.dirname+"[build_modules]: list_of_module_files in murat/"+self.subdir," : ",list_of_module_files

        for module in list_of_module_files:
            if (not module in skip_list):
                #    print ".cc file: "+cc
                o = '#/tmp/src/'+self.d1+'/'+string.split(module,'.')[0]+'.o'
                env.SharedObject(o,module)

                mname = string.split(os.path.basename(module),'.')[0];
                lib   = '#/lib/libmu2e_'+self.dirname+'_'+mname+'.so';
                #    print "o: "+o, "lib:"+lib
                env.SharedLibrary(lib,o,LIBS = ['libStntuple_mod.so',libs]);

# Export the class so that it can be used in the SConscript files
# For reasons even Rob does't understand, this must come before the env.SConscript(ss) line.
Export('murat_helper')

 
