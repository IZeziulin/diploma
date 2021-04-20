import par, kin
import numpy as np

global NsigmaX
NsigmaX = 20

def integrand1(x, p):
	"""p changes in range (0; 1) for incoherent case and (-inf; +inf) else"""
	x0 = par.mX - NsigmaX * par.GX
	xm = par.mX + NsigmaX * par.GX
	deltaX = xm - x0
	y0 = 2 * p / deltaX
	ym = 2 / deltaX - y0
	return (ym - y0) * (x - x0) / deltaX + y0

def integrand2(x, par0, par1):
	return np.sqrt(np.abs(kin.RBW(x, par0, par1))**2)
	
def integrandRBW(x, par0, par1):
	return np.abs(kin.RBW(x, par0, par1))**2
	
def integrand3(x, pars):
	from scipy.integrate import quad
	xMin = par.mX - NsigmaX * par.GX
	xMax = par.mX + NsigmaX * par.GX
	c2 = quad(integrand2, xMin, xMax, args = (pars[3], pars[4]))
	return np.abs(integrand1(x, pars[0]) + pars[1] * np.exp(1j*pars[2]) * kin.RBW(x, pars[3], pars[4]) / c2[0])**2
	
def fInCoher(x, pars):
	from scipy.integrate import quad
	import numpy as np
	xMin = par.mX - NsigmaX * par.GX
	xMax = par.mX + NsigmaX * par.GX
	c1 = quad(integrand1, xMin, xMax, args = (pars[0],))
	c2 = quad(integrandRBW, xMin, xMax, args = (pars[2], pars[3]))
	return pars[1] * integrand1(x, pars[0]) / c1[0] + (1 - pars[1]) * integrandRBW(x, pars[2], pars[3]) / c2[0]
    
def fCoher(x, pars):
	from scipy.integrate import quad
	xMin = par.mX - NsigmaX * par.GX
	xMax = par.mX + NsigmaX * par.GX
	c = quad(integrand3, xMin, xMax, args = pars)
	return integrand3(x, pars) / c[0]
	
def IntfInCoher(x, pars):
	from scipy.integrate import quad 
	xMin = par.mX - NsigmaX * par.GX
	xMax = par.mX + NsigmaX * par.GX
	c = quad(fInCoher, xMin, xMax, args = pars)
	return c[0]

def fitInCoher(data, NsigmaX):
	import iminuit
	
	xMin = par.mX - NsigmaX * par.GX
	xMax = par.mX + NsigmaX * par.GX
	init = [0.2, 0.5, 3872, 1]
	d_limits = [0, 0, xMin, 0]
	u_limits = [1, 1, xMax, xMax - xMin]

	def loglh(a1, a2, a3, a4):
		return - np.sum(np.log(fInCoher(data, [a1, a2, a3, a4])))
		
	minimizer = iminuit.Minuit(loglh, a1 = init[0], a2 = init[1], a3 = init[2], a4 = init[3])
	minimizer.errordef = iminuit.Minuit.LIKELIHOOD
	minimizer.fixed = (False, False, True, False)
	for i in range (0, 4):
		minimizer.limits[i] = [d_limits[i], u_limits[i]]
	minimizer.migrad()
	return list(minimizer.values)
	
def drawFit(data, f_name, f_val, N_bins, ax):
	import matplotlib.pyplot as plt
	x = np.linspace(xMin, xMax, N_bins)
    delta = (xMax - xMin) / N_bins
    ax.hist(data_all_err[i], bins = N_bins)
    nc = len(data_all_err[i]) * delta
    ax.plot(x, nc*fitter.fInCoher(x, fit_values[i]))	

