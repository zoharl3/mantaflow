
# manta probably overrides python's cout, and tee messes it

import os, sys

cmd = r'..\build\debug\manta zflip.py'
if 1: # release
    cmd = r'..\build\RelWithDebInfo\manta zflip.py'

cmd += r' 2>&1 | python \prj\python\tee.py _log.ans'

os.system( cmd )

os.system( "copy_log.bat" )

print( 'run.py is done' )
