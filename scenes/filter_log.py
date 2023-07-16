
import os, sys, re

sys.path.append( r'c:\prj\python\\' )
from text_color import *

os.system( 'cls' )

f = open( '_log.ans' )

lines = f.readlines()

n_print = 0 # number of lines to print without test
i = 0
while i < len( lines ):
    line = lines[i].rstrip()
    #print ( f'i={i}: { line }' )
    
    prn = 0
    if n_print > 0:
        prn = 1
        n_print -= 1
        
    elif re.search( '^- time:', line ):
        print()
        emphasize( line )
        
    elif re.search( '- push_particles', line ):
        prn = 1
        n_print = 1
        
    elif re.search( 'push time:', line ):
        prn = 1
        
    if prn:
        print( line )
        
    i += 1

input()
