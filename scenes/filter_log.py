
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
        'push time:',
        'gurobi optimize:'
        ]
    
    prn = 0
    if print_next > 0:
        prn = 1
        print_next -= 1
        
    elif re.search( '^- time:', line ):
        print()
        emphasize( line )
        
    elif re.search( '- push_particles', line ):
        prn = 1
        print_next = 1
        
    elif any( s in line for s in ss ):
        prn = 1
        
    if prn:
        print( line )
        
    i += 1

input()