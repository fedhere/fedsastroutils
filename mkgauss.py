import sys
import numpy as np
import optparse


def mkgaussian(size, center=None,sig=None,theta=None):
    x = np.arange(size*1.0)
    y = x[:,np.newaxis]
    if sig==None:
        sig=(size/10,size/10)
    elif isinstance( sig, ( int, long ) )  or isinstance( sig, ( float ) ):
        sig=(sig,sig)
    else:
        assert len(sig)==2 ,"sig must be an integer/float, a list/array of 2 (or None)"
    if center is None:
        x0 = y0 =0.5*size
    else :
        assert (len(center)==2) ,"center must be a list/array of 2 (or None)"
        x0 = center[0]
        y0 = center[1]

    if not theta:
        return np.exp(-0.5*(((x-x0)/sig[0])**2 + ((y-y0)/sig[1])**2)),[x0,y0],sig,theta
    else:
        assert (isinstance( theta, ( float ) ) or isinstance( theta, ( int, long ) )), "theta must be a float"
        thetarad=theta*np.pi/180.
        a =  0.5*(np.cos(thetarad)/sig[0])**2 + 0.5*(np.sin(thetarad)/sig[1])**2
        b = -0.25*np.sin(2*thetarad)/sig[0]**2 + 0.25*np.sin(2*thetarad)/sig[1]**2 
        c =  0.5*(np.sin(thetarad)/sig[0])**2 + 0.5*(np.cos(thetarad)/sig[1])**2
        
        return  np.exp( - (a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2)),[x0,y0],sig,theta




if __name__=='__main__':
    parser = optparse.OptionParser(usage="python mkgauss.py ", conflict_handler="resolve")
    parser.add_option('--size', default=100, type=int,
                      help='size of the array (side)')
    parser.add_option('--center', default=None, type=str,
                      help='center of the gaussian: float or comma separated sequence of 2 floats')
    parser.add_option('--sigma', default=None, type=str,
                      help='sigma of the gaussian: float or comma separated sequence of 2 floats')
    parser.add_option('--theta', default=None, type=float,
                      help='angle in degreed')
    parser.add_option('--normalization', default=1.0, type=float,
                      help='normalization factor')
    parser.add_option('--save', default=False,action='store_true',
                      help='save the gaussian as a fits file')
    parser.add_option('--show', default=False,action='store_true',
                      help='show the gaussian surface')
    sigma,center,theta=None,None,None
    options,  args = parser.parse_args()
    if options.sigma:
        try:
            sigma=float(options.sigma)
        except:
            try:
                sigma=[float(s) for s in options.sigma.split(',')]
                assert len(sigma)==2 ,'sigma must be float or comma separated sequence of 2 floats'
            except:
                print 'sigma must be float or comma separated sequence of 2 floats'
                sys.exit()
    if options.center:
        try:
            center=float(options.center)
        except:
            try:
                center=[float(s) for s in options.center.split(',')]
                assert len(center)==2, 'center must be float or comma separated sequence of 2 floats'
            except:
                print 'center must be float or comma separated sequence of 2 floats'
                sys.exit()
        
        

    gauss,c,s,t=mkgaussian(options.size,sig=sigma,center=center,theta=options.theta)
    gauss = options.normalization*gauss
    if options.show:
        import pylab as pl
        pl.imshow(gauss)
        pl.show()
    if options.save:
        import pyfits 
        if t == None: t = 0
        hdu = pyfits.PrimaryHDU(gauss)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto('gaussian_center%.1fx%.1f_sigma%.1fx%.1f_theta%.1f_size%d.fits'
                        %(c[0],c[1],s[0],s[1],t,options.size),clobber=True)

    
