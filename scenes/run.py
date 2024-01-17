
import os, sys, keyboard, subprocess, shutil
from pathlib import Path

sys.path.append( r'c:\prj\python\\' )
import FileLock

lock = FileLock.FileLock( '/tmp/manta_run.py.lock' )

#   0       1          2          3           4
# FLIP, FIXED_VOL, CORRECT19, DE_GOES22, MATLAB_FLIP
#
methods = [0]
#methods = [0,1]
#methods = [2,0]
#methods = [2,1]
#methods = [0,2,1]
#methods = [1,0,2]
#methods = [0,2,1,3]

exe = r'..\build\debug\manta' # cmd can't execute "../"
if 1: # release
    exe = r'..\build\RelWithDebInfo\manta'

script = 'zflip.py'
out_dir_root = r'c:/prj-external-libs/mantaflow/out/'

cmd_base = f'"{exe}" {script}'

def run( method ):
    cmd = f'{cmd_base} {method}'

    os.system( 'rm _log.ans' )

    # shell redirection loses ascii color?
    # ConEmu supports it, but it's slow
    # https://stackoverflow.com/questions/3515208/can-colorized-output-be-captured-via-shell-redirect
    # https://stackoverflow.com/questions/12573574/redirect-command-prompt-output-to-gui-and-keep-color
    #cmd += r'> _log.ans'
    #cmd += r' 2>&1 | python \prj\python\tee.py _log.ans'
    #cmd = f'script --flush --quiet --return _log.ans --command "{cmd}"' # if the console is terminated, the cmd still runs?
    cmd += r' 2>&1 | tee _log.ans'

    print( cmd )
    #subprocess.run( cmd ) # runs a process instead of executing a shell command; an issue for tee
    os.system( cmd )
    print()

def main():
    ret = 0
    
    # delete first level dirs
    for path in Path( out_dir_root ).glob('*'):
        if path.is_dir():
            #print( f'Deleting "{path.name}"' )
            for path2 in path.glob("*"):
                if path2.is_file():
                    path2.unlink()
                else:
                    print( f'Error "{path2.is_file()}" isn\'t a file' )
                    return -1
            try:
                path.rmdir()
            except:
                print( f"Can't delete {path}" )
                return -1
        if path.is_file():
            path.unlink()

    # run
    for method in methods:
        run( method )

        # check if there's a log in the latest dir, which means it ended gracefully
        dirs = [ f for f in Path( out_dir_root ).iterdir() if f.is_dir() ]
        if not dirs:
            print( 'Error: no directories' )
            os.system( f'copy_log.bat "{out_dir_root}"' )
            break
        latest_dir = max( dirs, key=os.path.getmtime )
        log = latest_dir / '_log.ans'
        if not log.exists():
            print( f"log doesn't exist (premature/forced/user exit): '{log}'" )
            os.system( f'copy_log.bat "{latest_dir.as_posix()}"' )
            
            path = Path( r'c:\prj\test_data\relative\_tmp\_matlab_out.txt' )
            if path.exists():
                print( 'moving matlab log' )
                shutil.move( path.as_posix(), latest_dir.as_posix() )
            
            #ret = -1
            break

    print( 'run.py is done' )
    return ret

# __main__
if __name__ == "__main__":
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
