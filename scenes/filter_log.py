
import os, sys, re

sys.path.append( r'c:\prj\python\\' )
from text_color import *

os.system( 'cls' )

f = open( '_log.ans' )

lines = f.readlines()

print_next = 0 # number of lines to print without test
i = 0
while i < len( lines ):
    line = lines[i].rstrip()
    #print ( f'i={i}: { line }' )
    
    ss = [
        #'vel.MaxAbs',
        'gurobi optimize:',
        #'q_stat',
        #'push time:',
        'avg_opt_time',
        ]
    
    prn = 0
    
    # print next
    if print_next > 0:
        prn = 1
        print_next -= 1
        
    # time
    elif re.search( '^- time:', line ):
        print()
        emphasize( line )
        
    elif 0 and re.search( '- push_particles', line ):
        prn = 1
        print_next = 1
        
    # ss
    elif any( s in line for s in ss ):
        prn = 1
        
    if prn:
        print( line )
        
    i += 1

print( '\nEnter...' )
input()
