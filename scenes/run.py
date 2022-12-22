
import os, sys, keyboard

cmd = r'..\build\debug\manta zflip.py'
if 1: # release
    cmd = r'../build/RelWithDebInfo/manta zflip.py'
    #cmd = r'..\build\RelWithDebInfo\manta zflip.py'

os.system( 'rm _log.ans' )

# shell redirection loses ascii color; ConEmu supports it, but it's slow
# https://stackoverflow.com/questions/3515208/can-colorized-output-be-captured-via-shell-redirect
# https://stackoverflow.com/questions/12573574/redirect-command-prompt-output-to-gui-and-keep-color
#cmd += r'> _log.ans'
#cmd += r' 2>&1 | python \prj\python\tee.py _log.ans'
cmd = f'script --flush --quiet --return _log.ans --command "{cmd}"'

os.system( cmd )

os.system( "copy_log.bat" )

print( 'run.py is done' )

# pause
if 1:
    print( 'press esc or enter...' )
    while 1:
        ch = keyboard.read_key()
        if ch == 'esc' or ch == 'enter':
            break
elif 0:
    print( 'press a key...' )
    keyboard.read_key()
elif 0:
    print( 'press enter...' )
    input()
    