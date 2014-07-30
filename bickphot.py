#!/usr/bin/env python

import sys, numpy, scipy.special as special
import pyfits 
def gaussianImage(n, center_x, center_y, sigma):
    
    image = numpy.zeros([n, n], dtype=float)
    
    for ix in range(n):
        
        for iy in range(n):
            
            x, y = ix - center_x, iy - center_y
            
            A = 1.0/(2.0*numpy.pi*sigma**2)
            
            image[iy,ix] = A*numpy.exp(-(x**2 + y**2)/(2.0*sigma**2))

    out_fits = pyfits.PrimaryHDU(image) 
    out_fits.writeto('fitsfilehere.fits', clobber=True)            
    return image

def wijCoefficients(n, x, y, rad1, rad2, theta, ellipticity):
    
    wid, xcen, ycen = n, (n+1)/2, (n+1)/2
    
    dx, dy = x - xcen, y - ycen
    
    if n%2:
        
        dx, dy = dx + 1, dy + 1
        
    ftWij = numpy.zeros([wid, wid], dtype=complex)

    scale = 1.0 - ellipticity
        
    for iy in range(wid):

        ky = float(iy - ycen)/wid

        for ix in range(wid):

            kx = float(ix - xcen)/wid
            
# rotate and rescale
            
            cosT, sinT = numpy.cos(theta), numpy.sin(theta)
            
            kxr, kyr = kx*cosT + ky*sinT, scale*(-kx*sinT + ky*cosT)
            
            k = numpy.sqrt(kxr**2 + kyr**2)

# compute the airy terms, and apply shift theorem

            if k != 0.0:

                airy1 = rad1*special.j1(2.0*numpy.pi*rad1*k)/k
                
                airy2 = rad2*special.j1(2.0*numpy.pi*rad2*k)/k
                
            else:

                airy1, airy2 = numpy.pi*rad1**2, numpy.pi*rad2**2
                
            airy = airy2 - airy1
            
            phase = numpy.exp(-1.0j*2.0*numpy.pi*(dx*kxr + dy*kyr))
            
    ftWij[iy,ix] = phase*scale*airy

    ftWijShift = numpy.fft.fftshift(ftWij)
    
    wijShift = numpy.fft.ifft2(ftWijShift)
    
    wij = numpy.fft.fftshift(wijShift)
    
    return wij.real

if __name__ == '__main__':

    import pylab as pl
    n = int(sys.argv[1])
    
    x, y, psf_sigma = map(float, sys.argv[2:])
    
    radius_inner, theta, ellipticity = 0.0, numpy.pi*(0.0)/180.0, 0.0
    
    psf = gaussianImage(n, x, y, psf_sigma)
    
    growthf=numpy.zeros((len(numpy.arange(radius_inner + 0.1, 6.0*psf_sigma, 0.1)),3),float)
    
    for i,radius in enumerate(numpy.arange(radius_inner + 0.1, 6.0*psf_sigma, 0.1)):
        
        wij = wijCoefficients(n, x, y, radius_inner, radius, theta, ellipticity)
        
        flux_measured = (psf*wij).sum()
        
        flux_analytic = (1.0 - numpy.exp(-radius**2/(2.0*psf_sigma**2)))
        
        print radius, flux_measured, flux_analytic
        growthf[i][0],growthf[i][1],growthf[i][2]=radius, flux_measured, flux_analytic

    

    growthf=growthf.transpose()
    
    pl.plot(growthf[0],growthf[1],'r-')
#    pl.plot(growthf[0],growthf[2],'b-')
    pl.show()
