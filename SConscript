#!/usr/bin/env python
#
import os, re, string, subprocess
Import('env')
from stntuple_helper import *
#------------------------------------------------------------------------------
# last two components of the path. Ex: /not/this/but/THIS/AND_THIS
#                                      "AND_THIS" is usually "src"
#------------------------------------------------------------------------------
x = subprocess.call('scripts/build_config',shell=True)

murat_env = env.Clone()

murat_env['CPPPATH' ].append('-I'+os.environ['BUILD_BASE']+'/include');
murat_env['CXXFLAGS'].append('-I'+os.environ['BUILD_BASE']+'/include');
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
murat_env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen})
murat_env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('murat_env')
Export('stntuple_helper')
