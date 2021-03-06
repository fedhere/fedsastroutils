import numpy as np

def bayesian_blocks(t, fp_rate=.05, BINNED=False):
    '''Bayesian Blocks Implementation based on Scargle's 2013 paper 
    http://adsabs.harvard.edu/abs/2012arXiv1207.5578S
    and  Jake Vanderplas python  implementation 
    https://jakevdp.github.io/blog/2012/09/12/dynamic-programming-in-python/

    but it is implemented for time series (the histogramming functionality 
    which is most useful for time tagged events etc 
    is maintaned through the option BINNED=True)
   
    t : if BINNED: 1darray,        data to be histogrammed
       else: 2darray, the first dimension is interpreted as 
                      timestamps of the measurementm which are 
                      assumed to be regularly spaced (for now)

    fp_rate: false positive rate, default 5%

    Returns
    -------
    bins : ndarray
        array containing the N bin edges where N and the edge locations 
        are optimal in a bayesian_block sense

    use in conjunction with pylabsetup.myhisteps to get steps to plot
    ''' 


    if BINNED:
        tt=np.sort(t)
        num_points = tt.size
        edges = np.concatenate([tt[:1],0.5*(tt[1:] + tt[:-1]), tt[-1:]])
        block_length = tt[-1] - edges        
        ncp_prior = 4.0-np.log(fp_rate/(0.0136*num_points**0.478))
        #print ncp_prior
    else:
        num_points = t[0].size
        tt=np.arange(num_points)#np.arange(num_points)
        tones=np.ones(num_points)
        edges = np.concatenate([tt[:1],0.5*(tt[1:] + tt[:-1]), tt[-1:]])
        block_length = tt[-1] - edges
        
        nn_vec=[]
        #block_length = t[0][-1] - edges
        #print fp_rate
        ncp_prior = fp_rate/0.05*(1.1 - 0.6 * np.log10(num_points))
#(1.32)/fp_rate + 0.577 * np.log10(num_points)
#4.0-np.log(fp_rate/(0.0136*num_points**0.478))
#fp_rate#np.log10(num_points)*0.6+1.3
        #

    nn_vec = np.ones(num_points)
    best = np.zeros(num_points, dtype=float)
    last = np.zeros(num_points, dtype=int)
    #-----------------------------------------------------------------
    # Start with first data cell; add one cell at each iteration
    #-----------------------------------------------------------------

    if BINNED:
        for K in range(num_points):
            # Compute the width and count of the final bin for all possible
            # locations of the K^th changepoint
            width = block_length[:K + 1] - block_length[K + 1]
            count_vec = np.cumsum(nn_vec[:K + 1][::-1])[::-1]

            # evaluate fitness function for these possibilities
            fit_vec = count_vec * (np.log(count_vec) - np.log(width))
            fit_vec -= ncp_prior  # 4 comes from the prior on the number of changepoints
            fit_vec[1:] += best[:K]

            # find the max of the fitness: this is the K^th changepoint
            i_max = np.argmax(fit_vec)
            last[K] = i_max
            best[K] = fit_vec[i_max]

    else:
        for K in range(num_points):
            #print K+1,
            sum_1=np.cumsum(t[1][:K + 1][::-1])
            sum_0=np.cumsum(tones[:K + 1][::-1])
            #print sum_0[:K + 1][::-1][0],sum_1[:K + 1][::-1][0],
            fit_vec=((sum_1[:K + 1][::-1])**2)/(4.0*sum_0[:K + 1][::-1])
            #print fit_vec[0]
            fit_vec -= ncp_prior
            fit_vec[1:] += best[:K]

            # find the max of the fitness: this is the K^th changepoint
            i_max = np.argmax(fit_vec)
            last[K] = i_max
            best[K] = fit_vec[i_max]
    

    change_points =  np.zeros(num_points, dtype=int)
    i_cp = num_points
    ind = num_points
    while True:
        i_cp -= 1
        change_points[i_cp] = ind
        if ind == 0:
            break
        ind = last[ind - 1]
    change_points = change_points[i_cp:]

    #print change_points
    if not BINNED: return edges[change_points]/(max(edges[change_points])-min(edges[change_points]))*(max(t[0])-min(t[0]))+min(t[0])
    #print "hallo"
    return edges[change_points]
