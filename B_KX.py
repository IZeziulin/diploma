import kin, par
import numpy as np

global Nsigma
Nsigma = 10

"""
Variables:
0. thetaX
1. thetaJPsi
2. ThetaPiPi
3. phiJPsi
4. phiPiPi
5. mPiPi
6. mX (mJPsiPiPi)
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
    pX = kin.p_decay12_d1(par.mB, data[6], par.mK)
    return 10**9 * pX * kin.RBW(data[6], par.mX, par.GX) * kin.FormFactorL1(pX) *\
    sgn * kin.wignerD(data[0], 0, 1, 0, lmbdJPsi-lmbdPiPi) * pPi * \
    kin.wignerD(data[1], data[3], 1, lmbdJPsi, lmbdDltMu) *\
    kin.FormFactorL1(pPi) * kin.wignerD(data[2], data[4], 1, lmbdPiPi, 0) *\
    kin.RBW(data[5], par.mRho, kin.GRhoFn(data[5])) 

def M_Xmod(data, lmbdJPsi, lmbdPiPi, lmbdDltMu):
    """returns matrix element for X->J/PsiPiPi"""
    if lmbdJPsi+lmbdPiPi == 0:
        return 0
    if lmbdJPsi+lmbdPiPi > 0:
        sgn = 1
    else:
        sgn = -1
    pPi = kin.p_decay12_d1(data[5], par.mPi, par.mPi)
    pX = kin.p_decay12_d1(par.mB, data[6], par.mK)
    return 10**9 * pX * kin.T_in(1, data[6] - par.m_00st) * kin.FormFactorL1(pX) *\
    sgn * kin.wignerD(data[0], 0, 1, 0, lmbdJPsi-lmbdPiPi) * pPi * \
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
                s += M_Xmod(data, lmbdJPsi, lmbdPiPi, lmbdDltMu)        
        t += abs(s)**2
    return t * kin.p_decay12_d1(par.mB, data[6], par.mK) * kin.p_decay12_d1(data[6], data[5], par.mJPsi) *\
           kin.p_decay12_d1(data[5], par.mPi, par.mPi) * np.sin(data[0]) * np.sin(data[1]) * np.sin(data[2]) * data[5] / (data[5] ** 2 * data[6] ** 2)

def GenerVariables(N):
    """random generation of variables for further calculations"""
    data = np.zeros(7 * N)
    data.shape = 7, N
    for i in range (0, 3):
        data[i] = np.pi * np.random.random(N)
    for i in range (3, 5):
        data[i] = 2 * np.pi * np.random.random(N)
    #Nsigma = 1.5
    data[6] = par.mX + Nsigma * par.GX * (-1 + 2 * np.random.random(N))
    data[5] = 2 * par.mPi + (par.mX + Nsigma * par.GX - par.mJPsi - 2 * par.mPi) * np.random.random(N)
    check = data[6] - par.mJPsi
    mask = data[5] < check
    dim = len(data[0][mask])
    datanew = np.zeros(7 * dim)
    datanew.shape = 7, dim
    for i in range (0, 7):
        datanew[i] = data[i][mask]
    return datanew

def GenerData(N):
    """process data generation"""
    variables = GenerVariables(N)
    data = G(variables)
    xi = data.max() * np.random.random(len(data))
    print(data.max())
    mask = data > 1.05 * xi
    dim = len(variables[0][mask])
    varnew = np.zeros(7 * dim)
    varnew.shape = 7, dim
    for i in range (0, 7):
        varnew[i] = variables[i][mask]
    return varnew

def test_thetaX_to_mKPiPi(data):
	pK = kin.p_decay12_d1(par.mB, par.mK, data[6])
	pJPsi_Xframe = kin.p_decay12_d1(data[6], par.mJPsi, data[5])
	beta = kin.Beta(data[6], pK)
	#pJPsix = kin.BoostPx(par.mJPsi, pJPsi_Xframe, data[0], beta)
	#pJPsi = np.sqrt(pJPsix**2 + pJPsi_Xframe**2 * np.sin(data[0])**2)
	#EJPsi = kin.E(par.mJPsi, pJPsi)
	#EKPiPi = par.mB - EJPsi
	#print(EJPsi)
	EJPsi = kin.BoostE(par.mJPsi, pJPsi_Xframe, data[0], beta)
	#print(EJPsi2)
	#mKPiPi = np.sqrt(EKPiPi**2 - pJPsi**2)
	#print(mKPiPi)
	mKPiPi = np.sqrt(par.mB**2 + par.mJPsi**2 - 2* par.mB * EJPsi)
	#print(mKPiPi2) 
	return mKPiPi	





