def pool(data, minval, maxval, func, Nbins):
    import numpy as np
    expected = np.zeros(Nbins)
    pool = np.zeros(Nbins)
    experimental, t = np.histogram(data, Nbins) 
    w = sum(experimental)
    delta = (maxval - minval) / Nbins
    for i in range(0, Nbins):
        expected[i] = func((i + 0.5) * delta) * w 
    pool = (experimental - expected)**2 / np.sqrt(experimental)
    return pool 
  
def test_chi2_pool(data, minval, maxval, func, Nbins, ax):    
    import numpy as np
    pool_pr = pool(data, minval, maxval, func, Nbins)
    chi2 = sum(pool_pr)
    ax.hist(data, Nbins, color='r')
    ax.hist(pool_pr*500, Nbins, color='b')
