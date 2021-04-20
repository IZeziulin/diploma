"""
Common variables:
0. mKPiPi
1. thetaJPsi
2. ThetaPiPi
3. phiJPsi
4. phiPiPi
5. mPiPi
6. mJPsiPiPi
"""

"""
Parameters:
"0" - amplitude, KX
"1" - phase, K1JPsi
"2" - amplitude, K1JPsi
"""

import kin, par
import numpy as np

global NsigmaX
NsigmaX = 10

def FormDataKX(data):
    """makes variables for decay B->XK"""
    newdata = np.zeros(7 * len(data[0]))
    newdata.shape = 7, len(data[0])
    for i in range (1, 7):
        newdata[i] = data[i]
    """thetaX calculation"""
    pK = kin.p_decay12_d1(par.mB, par.mK, data[6])
    mKJPsi = np.sqrt(par.mK**2 + par.mB**2 + data[5]**2 + par.mJPsi**2 - data[0]**2 - data[6]**2)
    pPiPi = kin.p_decay12_d1(par.mB, mKJPsi, data[5])
    pMuMu = kin.p_decay12_d1(par.mB, par.mJPsi, data[0])
    costheta1 = (pK**2 + pPiPi**2 - pMuMu**2) / (2 * pPiPi * pK)
    beta = kin.Beta(data[6], -pK)
    pnewPiPix = kin.BoostPx(data[5], pPiPi, np.arccos(costheta1), beta)
    newdata[0] = np.pi - np.arccos(pnewPiPix / np.sqrt(pnewPiPix**2 + pPiPi**2 * (1 - costheta1**2)))
    #print(pK)
    #print(pPiPi)
    #print(pMuMu)
    return newdata
    
def FormDataK1JPsi(data):
    """makes variables for decay B->JPsiK1"""
    newdata = np.zeros(7 * len(data[0]))
    newdata.shape = 7, len(data[0])
    for i in range (1, 6):
        newdata[i] = data[i]
    newdata[6] = data[0]    
    """thetaK1 calculation"""
    pK = kin.p_decay12_d1(par.mB, par.mK, data[6])
    mKJPsi = np.sqrt(par.mK**2 + par.mB**2 + data[5]**2 + par.mJPsi**2 - data[0]**2 - data[6]**2)
    pPiPi = kin.p_decay12_d1(par.mB, mKJPsi, data[5])
    pMuMu = kin.p_decay12_d1(par.mB, par.mJPsi, data[0])
    costheta1 = (pMuMu**2 + pPiPi**2 - pK**2) / (2 * pPiPi * pMuMu)
    beta = kin.Beta(data[0], -pMuMu)
    pnewPiPix = kin.BoostPx(data[5], pPiPi, np.arccos(costheta1), beta)
    newdata[0] = np.pi -  np.arccos(pnewPiPix / np.sqrt(pnewPiPix**2 + pPiPi**2 * (1 - costheta1**2)))
    return newdata

def G(data, param):
    """returns total width for interference"""
    import B_KX as gen1
    import B_K1JPsi as gen2
    dataKX = FormDataKX(data)
    dataK1JPsi = FormDataK1JPsi(data)
    t = 0
    for lmbdDltMu in [-1, 1]:
        s = 0
        for lmbdJPsi in [-1, 0, 1]:
            for lmbdPiPi in [-1, 0, 1]:
                s += param[0] * gen1.M(dataKX, lmbdJPsi, lmbdPiPi, lmbdDltMu) + gen2.M(dataK1JPsi, lmbdJPsi,\
                     lmbdPiPi, lmbdDltMu) * np.exp(1j * param[1]) * param[2]
        t += abs(s)**2 
    return t * kin.p_decay12_d1(data[5], par.mPi, par.mPi) * np.sin(data[1]) * np.sin(data[2]) * data[0] * data[6] * data[5]/ data[5]**2

def GenerVariables(N):
    """random generation of variables for further calculations"""
    import dalitz
    #NsigmaX = 10
    data = np.zeros(7 * N)
    data.shape = 7, N
    """generation mJPsiPiPi"""
    data[6] = par.mX - NsigmaX * par.GX + 2 * NsigmaX * par.GX * np.random.random(N) 
    """generation mPiPi"""
    data[5] = 2 * par.mPi + (par.mB - par.mK - par.mJPsi - 2 * par.mPi) * np.random.random(N)
    """generation mKPiPi"""
    data[0] = par.mK + 2 * par.mPi + (par.mB - par.mJPsi - par.mK - 2 * par.mPi) * np.random.random(N)
    for i in range (1, 3):
        data[i] = np.pi * np.random.random(N)
    for i in range (3, 5):
        data[i] = 2 * np.pi * np.random.random(N)
    """check if mJPsiPiPi satisfies"""
    mask1 = (data[5] + par.mJPsi < data[6]) & (data[6] < par.mB - par.mK)
    datanew1 = np.zeros(7 * len(data[0][mask1]))
    datanew1.shape = 7, len(data[0][mask1])
    for i in range (0, 7):
        datanew1[i] = data[i][mask1]
    #print(len(datanew1[0]))    
    """check if mKPiPi satisfies"""
    mask2 = (dalitz.s23min(par.mB, datanew1[6], par.mJPsi, datanew1[5], par.mK) < datanew1[0]**2) & \
            (datanew1[0]**2 < dalitz.s23max(par.mB, datanew1[6], par.mJPsi, datanew1[5], par.mK))
    datanew2 = np.zeros(7 * len(datanew1[0][mask2]))
    datanew2.shape = 7, len(datanew1[0][mask2])
    for i in range (0, 7):
        datanew2[i] = datanew1[i][mask2]
    #print(len(datanew2[0]))    
    return datanew2

def GenerData(variables, par):
    """process data generation"""
    import numpy as np
    data = G(variables, par)
    xi = data.max() * np.random.random(len(data))
    mask = data > xi
    dim = len(variables[0][mask])
    varnew = np.zeros(7 * dim)
    varnew.shape = 7, dim
    for i in range (0, 7):
        varnew[i] = variables[i][mask]
    return varnew
	
def MakeData(N:int, par:list, filename:str):
	f = open(filename, "w")
	counter = 0
	while counter < N:
		var = GenerVariables(10**7)
		data = GenerData(var, par)
		counter += len(data[0])
		print(counter, "completed", sep = " ")
		for i in range (0, len(data[0])):
			for j in range (0, 7):
				f.write(str(data[j][i]) + " ")
			f.write("\n")
	f.close()		
		
def ReadData(filename:str):
	f = open(filename, "r")
	data = [[],[],[],[],[],[],[]]
	while True:
		s = f.readline()
		if s == "":
			break
		else:
			s = s.rstrip("\n")
			x = s.split(' ')
			for i in range (0, 7):
				data[i].append(float(x[i]))
	return np.array(data)
		
def interference(data, param):
	import B_KX as gen1
	import B_K1JPsi as gen2
	dataKX = FormDataKX(data)
	dataK1JPsi = FormDataK1JPsi(data)
	t = 0
	for lmbdDltMu in [-1, 1]:
		s = 0
		s1 = 0
		s2 = 0
		for lmbdJPsi in [-1, 0, 1]:
			for lmbdPiPi in [-1, 0, 1]:
				s += param[0] * gen1.M(dataKX, lmbdJPsi, lmbdPiPi, lmbdDltMu) + gen2.M(dataK1JPsi, lmbdJPsi,\
					lmbdPiPi, lmbdDltMu) * np.exp(1j * param[1]) * param[2]
				s1 += param[0] * gen1.M(dataKX, lmbdJPsi, lmbdPiPi, lmbdDltMu)
				s2 += gen2.M(dataK1JPsi, lmbdJPsi, lmbdPiPi, lmbdDltMu) * np.exp(1j * param[1]) * param[2]
			t += abs(s)**2 - abs(s1)**2 - abs(s2)**2
	return t * kin.p_decay12_d1(data[5], par.mPi, par.mPi) * np.sin(data[1]) * np.sin(data[2]) * data[0] * data[6] * data[5] / data[5]**2









