import matplotlib.pyplot as plt
import numpy as np
import interference as genI

def makehist1D(fig, axes, nrows, ncols, ar_title, data, Nbins):
    for i in range (0, nrows):
        for j in range (0, ncols):
            axes[i, j].hist(data[i * ncols + j], bins = Nbins)
            axes[i, j].set(title = ar_title[i * ncols + j])

def drawFunc(ax, xMin, xMax, Nbins, len_array, func):
    x = np.linspace(xMin, xMax, Nbins)
    norm = len_array * (xMax - xMin) / Nbins
    y = norm * func(x)
    ax.plot(x, y, 'r')
    
def DrawHist1D(data, save = False, name = ''):
	Nbins = 100
	title = ['$\\cos\\theta_{X}$', '$\\cos\\theta_{J/\\psi}$', '$\\cos\\theta_{\\pi\\pi}$', '$\\phi_{J/\\psi}$', '$\\phi_{\\pi\\pi}$', 
	'$m_{\\pi\\pi}, GeV$', '$m_{J/\\psi\\pi\\pi}, GeV$', '$m_{K\\pi\\pi}, GeV$']
	fig, axes = plt.subplots(ncols = 3, nrows = 3, figsize = (24, 18))
	dataKX = genI.FormDataKX(data)
	for i in range (0, 3):
		dataKX[i] = np.cos(dataKX[i])
	for i in range (5, 7):
		dataKX[i] = dataKX[i] / 1000	
	for i in range (0, 7):
		axes[i//3][i%3].hist(dataKX[i], Nbins)
		axes[i//3][i%3].set_xlabel(title[i], fontsize = 18)
		axes[i//3][i%3].set_ylabel('Frequency', fontsize = 18)	  
	axes[2][1].hist(data[0] / 1000, Nbins)
	axes[2][1].set_xlabel(title[7], fontsize = 18)
	axes[2][1].set_ylabel('Frequency', fontsize = 18)  
	if save:
		fig.savefig(name)
	plt.show()
    
def Draw_mJPsi_on_all(data, save = False, name = ''):
	title = ['$\\cos\\theta_{X}$', '$\\cos\\theta_{J/\\psi}$', '$\\cos\\theta_{\\pi\\pi}$', '$\\phi_{J/\\psi}$', '$\\phi_{\\pi\\pi}$', 
	'$m_{\\pi\\pi}, GeV$', '$m_{J/\\psi\\pi\\pi}, GeV$']
	dataKX = genI.FormDataKX(data)
	for i in range (0, 3):
		dataKX[i] = np.cos(dataKX[i])
	for i in range (5, 7):
		dataKX[i] = dataKX[i] / 1000
	fig, axes = plt.subplots(nrows = 2, ncols = 3, figsize = (12, 6))   
	for i in range (0, 2):
		for j in range (0, 3): 
			axes[i][j].scatter(data[6], data[3*i + j], s = 0.1)
			axes[i][j].set_xlabel(title[6], fontsize = 18)
			axes[i][j].set_xlabel(title[3*i + j], fontsize = 18)
	if save:
		fig.savefig(name)
	 
def DrawAngles2D(data, save = False, name = ''):
	fig, ax = plt.subplots(nrows = 2, ncols = 5, figsize = (20 ,8))
	title = ['$\\cos\\theta_{X}$', '$\\cos\\theta_{J/\\psi}$', '$\\cos\\theta_{\\pi\\pi}$', '$\\phi_{J/\\psi}$', '$\\phi_{\\pi\\pi}$']
	dataKX = genI.FormDataKX(data)
	for i in range (0, 3):
		dataKX[i] = np.cos(dataKX[i])
	counter = 0
	for i in range (0, 5):
		for j in range (i , 5):
			counter += 1
			ax[counter//5][counter % 5].scatter(dataKX[i], dataKX[j], s = 0.1)
			ax[counter//5][counter % 5].set_xlabel(title[i], fontsize = 18)
			ax[counter//5][counter % 5].set_ylabel(title[j], fontsize = 18)
	if save:
		fig.savefig(name)
	plt.show()
