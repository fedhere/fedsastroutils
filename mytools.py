import numpy as np
import scipy.stats as stats

def maskarray(myarray, value):
    mymask=np.asarray([False]*len(myarray))
    mymask[np.where(myarray == value)]=True
    return mymask

def sigma_clip(x, sig=3, method='median', iterations=100, verbose=False):
    ''' sigma clipping  -  iterating
    over 1D data rejecting points more than 
    a specified number (sig) of
    standard deviations away from the mean/median.
    '''
    data = np.array(x, copy=False)
    data = data.ravel()

#    mask = np.ones(data.size, bool)
    for i in range(iterations):
        if method == 'median':
            tmpdata = data - stats.nanmedian(data)
        elif method == 'mean':
            tmpdata = data - stats.nanmean(data)            

        std=stats.nanstd(data)
        ldold=len(data)
        data=data[np.abs(tmpdata)<sig*std]
            
        if verbose:print data
        if ldold==len(data):
            return data
    print "warning: mmax iterations (%d) reached"%iterations
    return data
