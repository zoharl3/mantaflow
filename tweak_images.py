
# There's a break in the end of the loop to prevent it from recurse

import sys, os, re, shutil
import colorama

# I use backslash, since I get files with backslashes from the os
path = '.\\'

bRecursive = 0

colorama.init()

lf = []
for root, dirs, files in os.walk(path):
    for file in files:
        m = re.match('^(.+)\.(png|jpg|bmp)$', file)
        if m:
            ext = m.group(2)
            fl = os.path.join(root, file)
            print(colorama.Fore.GREEN + colorama.Style.BRIGHT + fl + colorama.Style.RESET_ALL)

            # replace color
            # --opaque, -fill: replace color of opaque with fill
            # -fuzz: acceptable distance from opaque color
            if 0:
                cmd = 'convert "%s" -fuzz 5%% -opaque "rgb(170,170,170)" -fill white "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
            
            # convert bmp to png (png is lossless, jpg isn't unless it's .jp2 -- jpg2000)
            if 0 and ext.lower() == 'bmp':
                dest = os.path.join(root, m.group(1) + '.png')
                print(colorama.Fore.GREEN + colorama.Style.BRIGHT + fl + " => " + dest + colorama.Style.RESET_ALL)
                cmd = "convert %s %s" % (fl, dest)
                os.system(cmd)
                fl = dest
            
            # convert images to low res .pdf
            if 0:
                dest = os.path.join(root, m.group(1) + '.pdf')
                print(colorama.Fore.GREEN + colorama.Style.BRIGHT + fl + " => " + dest + colorama.Style.RESET_ALL)
                #cmd = "convert %s %s" % (fl, dest)
                cmd = "convert %s -resample 16 %s" % (fl, dest)
                os.system(cmd)

            # convert to black and white
            if 0:
                cmd = 'convert "%s" -threshold 70%% "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
                
            # fill background in white
            if 0:
                # fuzz was 20
                cmd = 'convert "%s" -fuzz 2%% -fill white -draw "color 0,0 floodfill" "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
                
            # trim borders
            if 1:
                cmd = 'convert "%s" -trim "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
            
            # add white border on left and right
            if 0:
                cmd = 'convert "%s" -bordercolor white -border 10x0 "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)

            # add white border on all sides
            if 0:
                cmd = 'convert "%s" -bordercolor white -border 10 "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
                
            # set transparent background
            if 0:
                cmd = 'convert "%s" -transparent white "%s"' % (fl,fl)
                os.system(cmd)
                
            # density
            if 0:
                cmd = 'convert "%s" -units PixelsPerInch -density 72 "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
                
            # resample to low resolution (e.g. to create draft for latex)
            if 0:
                cmd = 'convert "%s" -resample 8 "%s"' % (fl,fl)
                #print cmd
                os.system(cmd)
                
    if not bRecursive:
        break   #prevent descending into subfolders

if 0:
    print('\nPress Enter...')
    input()
    
    