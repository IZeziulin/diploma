def makehist1D(fig, axes, nrows, ncols, ar_title, data, Nbins):
    for i in range (0, nrows):
        for j in range (0, ncols):
            axes[i, j].hist(data[i * ncols + j], bins = Nbins)
            axes[i, j].set(title = ar_title[i * ncols + j])

def drawFunc(ax, xMin, xMax, Nbins, len_array, func):
    import numpy as np
    x = np.linspace(xMin, xMax, Nbins)
    norm = len_array * (xMax - xMin) / Nbins
    y = norm * func(x)
    ax.plot(x, y, 'r')
