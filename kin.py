import par
import numpy as np

"""Functions for 1->2 decay"""

def Kallen(x, y, z):
    return (x - y - z)**2 - 4*y*z

def p_decay12_d1(M, m1, m2):
    """returns momentum of the daughter particle 1"""
    import numpy as np
    return np.sqrt(Kallen(M**2, m1**2, m2**2)) / (2 * M)

def E_decay12_d1(M, m1, m2):
    """returns energy of the faughter particle 1"""
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
    return Gamma(beta) * (p + beta * E(m,p) * np.cos(theta))

def BoostE(m, p, theta, beta):
    """from daughter to mother particle rest system"""
    return Gamma(beta) * (E(m,p) + beta * p * np.cos(theta))

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



