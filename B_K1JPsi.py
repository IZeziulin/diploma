import kin, par

"""
Variables:
0. thetaK1
1. thetaJPsi
2. ThetaPiPi
3. phiJPsi
4. phiPiPi
5. mPiPi
6. mK1 (mKPiPi)
"""

global Nsigma
Nsigma = 10

def M(data, lmbdJPsi, lmbdPiPi, lmbdDltMu):
    """returns matrix element for X->J/PsiPiPi"""
    if lmbdJPsi != 0:
        sgn = 1
    else:
        sgn = -1
    pPi = kin.p_decay12_d1(data[5], par.mPi, par.mPi)
    return 10**9 * sgn * kin.wignerD(data[0], 0, 1, lmbdJPsi, lmbdPiPi) * kin.RBW(data[6], par.mK1, par.GK1)  *\
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
                s += M(data, lmbdJPsi, lmbdPiPi, lmbdDltMu)        
        t += abs(s)**2
    return t * kin.p_decay12_d1(par.mB, data[6], par.mJPsi) * kin.p_decay12_d1(data[6], data[5], par.mK) *\
           kin.p_decay12_d1(data[5], par.mPi, par.mPi) * np.sin(data[0]) * np.sin(data[1]) * np.sin(data[2]) / (data[6] * data[5])**2 

def GenerVariables(N):
    """random generation of variables for further calculations"""
    import numpy as np
    data = np.zeros(7 * N)
    data.shape = 7, N
    for i in range (0, 3):
        data[i] = np.pi * np.random.random(N)
    for i in range (3, 5):
        data[i] = 2 * np.pi * np.random.random(N)
    #Nsigma = 1
    data[6] = par.mK1 + Nsigma * par.GK1 * (-1 + 2 * np.random.random(N))
    """check if mKPiPi satisfies (use 1->2)"""
    mask1 = data[6] < par.mB - par.mJPsi
    #data[5] = 2 * par.mPi + (par.mK1 + Nsigma * par.GK1 - par.mK - 2 * par.mPi) * np.random.random(N)
    data[5] = 2 * par.mPi + (par.mB - par.mK - par.mJPsi - 2 * par.mPi) * np.random.random(N)
    """check if mPiPi satisfies (use 1->2)"""
    mask2 = data[5] < data[6] - par.mK
    mask = mask1 * mask2
    dim = len(data[0][mask])
    datanew = np.zeros(7 * dim)
    datanew.shape = 7, dim
    for i in range (0, 7):
        datanew[i] = data[i][mask]
    return datanew

def GenerData(N):
    """process data generation"""
    import numpy as np
    variables = GenerVariables(N)
    data = G(variables)
    print(data.max())
    xi = data.max() * np.random.random(len(data))
    mask = data > 1.05 * xi
    dim = len(variables[0][mask])
    varnew = np.zeros(7 * dim)
    varnew.shape = 7, dim
    for i in range (0, 7):
        varnew[i] = variables[i][mask]
    return varnew
    
