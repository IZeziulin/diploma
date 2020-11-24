import kin, par

"""
Variables:
0. thetaX
1. thetaJPsi
2. ThetaPiPi
3. phiJPsi
4. phiPiPi
5. mPiPi
"""

def M(data, lmbdJPsi, lmbdPiPi, lmbdDltMu):
    """returns matrix element for X->J/PsiPiPi"""
    if lmbdJPsi+lmbdPiPi == 0:
        return 0
    if lmbdJPsi+lmbdPiPi > 0:
        sgn = 1
    else:
        sgn = -1
    pPi = kin.p_decay12_d1(data[5], par.mPi, par.mPi)
    return sgn * kin.wignerD(data[0], 0, 1, 0, lmbdJPsi-lmbdPiPi) * pPi * \
    kin.wignerD(data[1], data[3], 1, lmbdJPsi, lmbdDltMu) *\
    kin.FormFactorL1(pPi) * kin.wignerD(data[2], data[4], 1, lmbdPiPi, 0) *\
    kin.RBW(data[5], par.mRho, kin.GRhoFn(data[5])) 

def G(data):
    """returns total width for X->J/PsiPiPi"""
    import numpy as np
    t = 0
    for lmbdDltMu in [-1, 1]:
        s = 0
        for lmbdJPsi in [-1, 0, 1]:
            for lmbdPiPi in [-1, 0, 1]:
                s = s + M(data, lmbdJPsi, lmbdPiPi, lmbdDltMu)        
        t =  (t + abs(s)**2) 
    return t * kin.p_decay12_d1(data[5], par.mPi, par.mPi) * kin.p_decay12_d1(par.mX, data[5], par.mJPsi) *\
            np.sin(data[0]) * np.sin(data[1]) * np.sin(data[2])

def GenerVariables(N):
    """random generation of variables for further calculations"""
    import numpy as np
    data = np.zeros(6 * N)
    data.shape = 6, N
    for i in range (0, 3):
        data[i] = np.pi * np.random.random(N)
    for i in range (3, 5):
        data[i] = 2 * np.pi * np.random.random(N)
    data[5] = 2 * par.mPi + (par.mX - par.mJPsi - 2 * par.mPi) * np.random.random(N)
    return data

def GenerData(variables, N):
    """process data generation"""
    import numpy as np
    data = G(variables)
    xi = data.max() * np.random.random(N)
    mask = data > 1.05 * xi
    dim = len(variables[0][mask])
    varnew = np.zeros(6 * dim)
    varnew.shape = 6, dim
    for i in range (0, 6):
        varnew[i] = variables[i][mask]
    return varnew

