import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from myastrotools import offsets
import sys

ra_g,dec_g=float(sys.argv[1]),float(sys.argv[2])
ra1,dec1=float(sys.argv[3]),float(sys.argv[4])
print "Galaxy coords: ",ra_g, dec_g
print "Star coords: ",ra1, dec1

c1 = SkyCoord(ra_g*u.degree, dec_g*u.degree, frame='fk5')
c2 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='fk5')
pa=c1.position_angle(c2)
sep = c1.separation(c2)

print "separation (arcsec)",sep.arcsec, "PA",pa
off=offsets([[ ra_g,dec_g ]],[[ ra1,dec1 ]])
print "RA Dec offset (arcsec)", off[0]*360.,off[1]*360.
