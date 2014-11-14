import sys,os
import numpy as np

if len(sys.argv)<2:
    print "only argument: sextractor file name"
    sys.exit()
inname=sys.argv[1]
outname=open(inname+'.reg','w')

headercount,xcol,ycol,flux=0,None,None,None
f=open (inname,'r')
for l in f:
    if l.startswith('#'):
        if 'X_IMAGE' in l:
            xcol=int(l.split()[1])
        if 'Y_IMAGE' in l:
            ycol=int(l.split()[1])
        if 'FLUX_AUTO' in l:
            flux=int(l.split()[1])
        headercount+=1
    else : break

print xcol,ycol,headercount,flux
outname.write('''# Region file format: DS9 version 4.1                                            
                                               
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1                                                
                                                          
image
''')

if flux:
    cat=np.loadtxt(inname,skiprows=headercount, usecols=(xcol-1,ycol-1, flux-1))
    print cat.T
    norm=np.max(cat.T[2])/60.0

    for star in cat:
        print >> outname, "circle(%.3f,%.3f,%.3f)"%(star[0],star[1],np.max([star[2]/norm,5]))

else:
    cat=np.loadtxt(inname,skiprows=headercount, usecols=(xcol-1,ycol-1))
    for star in cat:
        print >> outname,"circle(%.3f,%.3f,%.3f)"%(star[0],star[1],20)




