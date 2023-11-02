
import os, sys, re
import matlab.engine

sys.path.append( r'c:\prj\python\\' )
from text_color import *

os.system( 'cls' )

f = open( '_log.ans' )

lines = f.readlines()

print_next = 0 # number of lines to print without test

# search strings
ss = [
    #'vel.MaxAbs',
    'gurobi optimize:',
    #'q_stat',
    #'push time:',
    'avg_opt_time',
    ]

# data vars to send to matlab
vtime = []
vss = [ [] for _ in range( len(ss) ) ]

t = -1
for i, line in enumerate( lines ):
    line = lines[i].rstrip()
    #print ( f'i={i}: { line }' )
    
    prn = 0
    
    # print next
    if print_next > 0:
        prn = 1
        print_next -= 1
        
    else:
        # time
        m = re.search( '^- time: (\d+(\.\d+)?)', line )
        if m:
            print()
            emphasize( line )
            
            t = float( m.group(1) )
            
            # save data
            if t > 0:
                vtime.append( t )
        
        # push
        elif 0 and re.search( '- push_particles', line ):
            prn = 1
            print_next = 1
            
        # ss
        else:
            for j, s in enumerate( ss ):
                if s in line:
                    prn = 1
                    
                    # data
                    if t > 0:
                        m = re.search( '(\d+(\.\d+)?)', line )
                        if m:
                            #print( m.group() )
                            x = float( m.group(1) )
                            #print( f'appending {x} to vss[{j}]' )
                            vss[j].append( x )
                            
                    break
        
    # print
    if prn:
        print( line )
        

# send vars
#print( vtime )
#print( len( vss[0] ), vss[0] )
if 1:
    #names = matlab.engine.find_matlab()
    #print( 'running matlab sessions:', names )
    
    #eng = matlab.engine.connect_matlab() # doesn't connect to a c++'s session
    eng = matlab.engine.start_matlab() # exiting the script kills the engine
    eng.desktop( nargout=0 )
    
    eng.workspace['vtime'] = vtime
    eng.workspace['vss'] = vss

# pause
print( '\nEnter...' )
input()
