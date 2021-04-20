import kin, par
import numpy as np

"""
M - mass of mother particle,
m1, m2, m3 - masses of daughter particles
"""

def s12(M, m3, E3):
    """s12 in M particle rest system"""
    import numpy as np
    return M**2 + m3**2 - 2 * M * E3

def E2rest(m12, m1, m2):
    return (m12**2 - m1**2 + m2**2) / (2 * m12)

def E3rest(M, m12, m3):
    return (M**2 - m12**2 - m3**2) / (2 * m12)

def s23min(M, m12, m1, m2, m3):
    """PDG, 48.23b"""
    return (E2rest(m12, m1, m2) + E3rest(M, m12, m3))**2 \
            - (np.sqrt(E2rest(m12, m1, m2)**2 - m2**2) + np.sqrt(E3rest(M, m12, m3)**2 - m3**2))**2

def s23max(M, m12, m1, m2, m3):
    """PDG, 48.23a"""
    return (E2rest(m12, m1, m2) + E3rest(M, m12, m3))**2 \
            - (np.sqrt(E2rest(m12, m1, m2)**2 - m2**2) - np.sqrt(E3rest(M, m12, m3)**2 - m3**2))**2

"""Drawing functions"""

def DrawDalitzScatter(m12, m23, M, m1, m2, m3, ax):
    ax.scatter(m12**2, m23**2, s = 0.1)
    x = np.linspace(m1 + m2, M - m3, len(m12))
    y1 = s23min(M, m12, m1, m2, m3)
    y2 = s23max(M, m12, m1, m2, m3)
    ax.plot(x**2, y1, 'r')
    ax.plot(x**2, y2, 'r')
    
def DrawDalitzScatterAll(data, save, name):
	"""Drawing all dalitz scatter plots in """
	fig, ax = plt.subplots(nrows = 1, ncols = 3, figsize = (22, 6))
	mKJPsi = np.sqrt(par.mB**2 + par.mK**2 + par.mJPsi**2 + data[5]**2 - data[0]**2 - data[6]**2)
	t_JPsiPiPi = '$m^{2}_{J/\\psi\\pi\\pi}, GeV^{2}$'
	t_KPiPi = '$m^{2}_{K\\pi\\pi}, GeV^{2}$'
	t_KJPsi = '$m^{2}_{KJ/\\psi}, GeV^{2}$'
	ax[0].scatter(data[6]**2 / 10**6, data[0]**2 / 10**6, s = 0.1)
	ax[0].set_xlabel(t_JPsiPiPi, fontsize = 18)
	ax[0].set_ylabel(t_KPiPi, fontsize = 18)
	ax[1].scatter(data[6]**2 / 10**6, mKJPsi**2 / 10**6, s = 0.1)
	ax[1].set_xlabel(t_JPsiPiPi, fontsize = 18)
	ax[1].set_ylabel(t_KJPsi, fontsize = 18)
	ax[2].scatter(data[0]**2 / 10**6, mKJPsi**2 / 10**6, s = 0.1)
	ax[2].set_xlabel(t_KPiPi, fontsize = 18)
	ax[2].set_ylabel(t_KJPsi, fontsize = 18)
	if save:
		fig.savefig(name)
	plt.show()



