import numpy as np
from scipy import interpolate

def cosmo_dist(z1, z2, h=None, omegaM=None, omegaL=None, lum=None, ang=None, \
                   Mpc=None):
    """
    NAME:
      cosmo_dist

    PURPOSE:
      Return a cosmological distance (comoving [default], luminosity,
      or angular diameter in m [default] or Mpc) between redshift z1
      and z2 for given omegaM and omegaL.

    CALLING SEQUENCE:
      dist = cosmo_dist(z1, z2, h=, omegaM=, omegaL=, lum=, ang=, Mpc= )

    INPUTS:
      z1 - lower redshift
      z2 - higher redshift

    OPTIONAL INPUTS:
      h      - dimensionless hubble parameter (default 0.7)
      omegaM - matter density of the Universe (default 0.3)
      omegaL - cosmological constant (default 0.7)

    KEYWORDS:
      lum - flag to return luminosity distance
      ang - flag to return angular diameter distance
      Mpc - flag to return distance in Mpc

    OUTPUTS:
      dist - cosmological distance

    OPTIONAL OUTPUTS:

    EXAMPLES:
      Run the included hoggxample() function to reproduce several
      plots from arXiv:astro-ph/9905116.

    COMMENTS:

    REVISION HISTORY:
      2013/03/11 - Written by Greg Dobler (KITP/UCSB)

    ------------------------------------------------------------
    """

# -------- Set defaults
    if h==None:      h      = 0.7 
    if omegaM==None: omegaM = 0.3
    if omegaL==None: omegaL = 0.7



# -------- Utilties
    omegak = 1.0 - omegaM - omegaL
    mpMpc  = 3.08568025e22 # [m/Mpc]
    H0     = 100.0*h # [km s^-1 Mpc^-1]
    sol    = 299792.458 # [km s^-1]
    DH     = sol/H0 # [Mpc]
    z2     = np.array(z2)
    nz2    = z2.size
    z      = z1 + (z2[nz2-1]-z1+0.1)*np.arange(0,1,1e-4) # 0.1 is for interpol
    dz     = z[1]-z[0]
    zp1    = z+1
    Ez     = np.sqrt(omegaM*zp1**3 + omegak*zp1**2 + omegaL)



# -------- Integrate the appropriate distance and interpolate (allows
#          arrays of z2's to be input)
    Dc = DH*dz*np.cumsum(1.0/Ez) # [Mpc]
    if z2<z.all():
        print z2,z,"not right"
        return [-1]
    else:
        Dc = interpolate.interp1d(z,Dc)
        Dc = Dc(z2)



# -------- Take into account non-zero curvature
    if abs(omegak) < 1e-4: 
        DM = Dc # [Mpc]
    elif omegak > 1e-4:
        DM = DH/np.sqrt(omegak) * np.sinh(np.sqrt(omegak)*Dc/DH) # [Mpc]
    elif omegak < -1e-4:
        DM = DH/np.sqrt(abs(omegak))*np.sin(np.sqrt(abs(omegak))*Dc/DH) # [Mpc]



# -------- Return the appropriate distance measure
    if lum!=None:
#        print 'COSMO_DIST: Returing luminosity distance...'
        dist = (1+z2)*DM
    elif ang!=None:
#        print 'COSMO_DIST: Returing angular diameter distance...'
        dist = 1/(1+z2)*DM
    else:
#        print 'COSMO_DIST: Returing comoving distance...'
        dist = DM



# -------- Convert to m and return
    return dist*mpMpc if Mpc==None else dist


#
# REPRODUCE A FEW PLOTS IN arXiv:astro-ph/9905116
#
def hoggxample():

# -------- utilities
    import matplotlib.pyplot as plt

    h   = 1.0
    H0  = 100.0*h # [km s^-1 Mpc^-1]
    sol = 299792.458 # [km s^-1]
    DH  = sol/H0 # [Mpc]
    z   = np.arange(0,5,0.05)


    # Figure 1
    #
    dist_eds = cosmo_dist(0.0,z,h=h,omegaM=1.0,omegaL=0.0,Mpc=True)
    dist_ld  = cosmo_dist(0.0,z,h=h,omegaM=0.05,omegaL=0.0,Mpc=True)
    dist_hl  = cosmo_dist(0.0,z,h=h,omegaM=0.2,omegaL=0.8,Mpc=True)

    plt.figure()
    plt.plot(z,dist_eds/DH,'-')
    plt.plot(z,dist_ld/DH,':')
    plt.plot(z,dist_hl/DH,'-.')
    plt.xlim([0,5])
    plt.ylim([0,3])
    plt.xlabel('redshift z', fontsize=15)
    plt.ylabel(r'proper motion distance $D_M/D_H$', fontsize=15)
    plt.show()


    # Figure 2
    #
    dist_eds = cosmo_dist(0.0,z,h=h,omegaM=1.0,omegaL=0.0,Mpc=True,ang=True)
    dist_ld  = cosmo_dist(0.0,z,h=h,omegaM=0.05,omegaL=0.0,Mpc=True,ang=True)
    dist_hl  = cosmo_dist(0.0,z,h=h,omegaM=0.2,omegaL=0.8,Mpc=True,ang=True)

    plt.figure()
    plt.plot(z,dist_eds/DH,'-')
    plt.plot(z,dist_ld/DH,':')
    plt.plot(z,dist_hl/DH,'-.')
    plt.xlim([0,5])
    plt.ylim([0,0.5])
    plt.xlabel('redshift z', fontsize=15)
    plt.ylabel(r'angular diameter distance $D_A/D_H$', fontsize=15)
    plt.show()


    # Figure 3
    #
    dist_eds = cosmo_dist(0.0,z,h=h,omegaM=1.0,omegaL=0.0,Mpc=True,lum=True)
    dist_ld  = cosmo_dist(0.0,z,h=h,omegaM=0.05,omegaL=0.0,Mpc=True,lum=True)
    dist_hl  = cosmo_dist(0.0,z,h=h,omegaM=0.2,omegaL=0.8,Mpc=True,lum=True)

    plt.figure()
    plt.plot(z,dist_eds/DH,'-')
    plt.plot(z,dist_ld/DH,':')
    plt.plot(z,dist_hl/DH,'-.')
    plt.xlim([0,5])
    plt.ylim([0,16])
    plt.xlabel('redshift z', fontsize=15)
    plt.ylabel(r'luminosity distance $D_L/D_H$', fontsize=15)
    plt.show()


    # Figure 4
    #
    dist_eds = cosmo_dist(0.0,z,h=h,omegaM=1.0,omegaL=0.0,Mpc=True,lum=True)
    dist_ld  = cosmo_dist(0.0,z,h=h,omegaM=0.05,omegaL=0.0,Mpc=True,lum=True)
    dist_hl  = cosmo_dist(0.0,z,h=h,omegaM=0.2,omegaL=0.8,Mpc=True,lum=True)

    plt.figure()
    plt.plot(z,5*np.log10(dist_eds/1e-5 * h),'-')
    plt.plot(z,5*np.log10(dist_ld/1e-5 * h),':')
    plt.plot(z,5*np.log10(dist_hl/1e-5 * h),'-.')
    plt.xlim([0,5])
    plt.ylim([40,50])
    plt.xlabel('redshift z', fontsize=15)
    plt.ylabel(r'distance modulus $D_M + 5 log h$ (mag)$', fontsize=15)
    plt.show()

    return
