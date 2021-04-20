import par
import numpy as np

"""Functions for 1->2 decay"""

def Kallen(x, y, z):
    return (x - y - z)**2 - 4*y*z

def p_decay12_d1(M, m1, m2):
    """returns momentum of the daughter particle 1"""
    return np.sqrt(Kallen(M**2, m1**2, m2**2)) / (2 * M)

def E_decay12_d1(M, m1, m2):
    """returns energy of the daughter particle 1"""
    return (M**2 - m2**2 + m1**2) / (2 * M)

"""Lorentz-transformation"""

def E(m, p):
    import numpy as np
    return np.sqrt(m**2 + p**2)

def Beta(m, p):
    import numpy as np
    return p / np.sqrt(m**2 + p**2)

def Gamma(beta):
    import numpy as np
    return 1 / np.sqrt(1 - beta**2)

def BoostPx(m, p, theta, beta): 
    """from daughter to mother particle rest system"""
    px = p * np.cos(theta)
    return Gamma(beta) * (px + beta * E(m,p))

def BoostE(m, p, theta, beta):
    """from daughter to mother particle rest system"""
    px = p * np.cos(theta)
    return Gamma(beta) * (E(m,p) + beta * px)

"""Other functions"""

def FormFactorL1(p):
    return np.sqrt(par.r**2 / (1 + (p*par.r)**2))

def RBW(m1, m2, G):
    return 1 / (m1**2 - m2**2 - 1j*m2*G)

def GRhoFn(mPiPi):
    pRho = p_decay12_d1(par.mRho, par.mPi, par.mPi)
    pPiPi = p_decay12_d1(mPiPi, par.mPi, par.mPi)
    return par.GRho * ( (pPiPi/pRho)**3 ) * (par.mRho/mPiPi) *\
    (FormFactorL1(pPiPi) / FormFactorL1(pRho))**2

def wignerD(theta, phi, j, m1, m2):
    from sympy.abc import x, y
    from sympy.utilities.lambdify import lambdify
    from sympy.physics.quantum.spin import Rotation as Wigner
    d = Wigner.D(j, m1, m2, y, x, -y).doit().evalf() 
    return lambdify([x, y], d)(theta, phi)

"""Functions for description X(3872) transition"""

def mu(m1, m2):
	return m1*m2/(m1+m2)
	
mu1 = mu(par.m_D0, par.m_Dst0)
mu2 = mu(par.m_Dp, par.m_Dstm)
	
def k(mu, E, delta):
	return np.sqrt((2+0*1j)*mu*(E - delta))	
	
def D(E):
	k1 = k(mu1, E, par.d1)
	k2 = k(mu2, E, par.d2)
	return par.g_s*par.g_t - k1*k2 + 1j*(par.g_s+par.g_t)*(k1+k2)/2	

def T_00(E):
	return ((par.g_s+par.g_t)/2 + 1j * k(mu1, E, par.d1)) / D(E)

def T_01(E):
	return (par.g_s - par.g_t) / (2*D(E))

def T_10(E):
	return (par.g_s - par.g_t) / (2*D(E))
	
def T_11(E):
	return ((par.g_s+par.g_t)/2 + 1j * k(mu2, E, par.d2)) / D(E)

def T_in(g, E):
	return T_00(E) + T_01(E) - g*(T_10(E) + T_11(E)) 
	
	
