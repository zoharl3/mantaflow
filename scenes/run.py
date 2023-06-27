
import os, sys, keyboard, subprocess
from pathlib import Path

methods = [1]
#methods = [0,1]
#methods = [0,2,1]

exe = r'../build/debug/manta'
if 1: # release
    exe = r'../build/RelWithDebInfo/manta'

script = 'zflip.py'
out_dir_root = r'c:/prj-external-libs/mantaflow/out/'

cmd_base = f'{exe} {script}'

def run( method ):
    cmd = f'{cmd_base} {method}'

    os.system( 'rm _log.ans' )

    # shell redirection loses ascii color?
    # ConEmu supports it, but it's slow
    # https://stackoverflow.com/questions/3515208/can-colorized-output-be-captured-via-shell-redirect
    # https://stackoverflow.com/questions/12573574/redirect-command-prompt-output-to-gui-and-keep-color
    #cmd += r'> _log.ans'
    #cmd += r' 2>&1 | python \prj\python\tee.py _log.ans'
    cmd = f'script --flush --quiet --return _log.ans --command "{cmd}"'

    result = subprocess.run( cmd )
    print( f'returncode={result.returncode}' )

def main():
    # delete first level dirs
    for path in Path( out_dir_root ).glob("*"):
        if path.is_dir():
            #print( f'Deleting "{path.name}"' )
            for path2 in path.glob("*"):
                if path2.is_file():
                    path2.unlink()
                else:
                    print( f'Error "{path2.is_file()}" isn\'t a file' )
                    return -1
            path.rmdir()

    # run
    for method in methods:
        run( method )

    print( 'run.py is done' )
    return 0

ret = main()

# pause
if 0 or ret != 0:
    if 1:
        print( '\nPress esc, space, or enter...' )
        while 1:
            ch = keyboard.read_key()
            if ch == 'esc' or ch == 'enter' or ch == 'space':
                break
    elif 0:
        print( '\nPress a key...' )
        keyboard.read_key()
    elif 0:
        print( '\nPress enter...' )
        input()
