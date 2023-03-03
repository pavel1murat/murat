#!/usr/bin/env python
#
import os, re, string, sys, subprocess
Import('env')
sys.path.append(os.getenv("MUSE_WORK_DIR")+'/site_scons')
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
# print("murat/SConscript: PWD:"+os.getenv("PWD"))
x = subprocess.call(os.getenv("MUSE_WORK_DIR")+'/murat/scripts/build_config_muse murat',shell=True)
# print("murat/SConscript: BACK")

murat_env = env.Clone()

murat_env['CPPPATH' ].append(os.environ['MUSE_WORK_DIR']+'/include');
murat_env.Append(FORTRANPATH = [os.environ['MUSE_WORK_DIR']+'/include']);
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
exec(open(os.environ['MUSE_WORK_DIR']+"/site_scons/stntuple_site_init.py").read())

from stntuple_helper import *

murat_env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen})
murat_env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('murat_env')
Export('stntuple_helper')
