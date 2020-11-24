import numpy as np

"""1d functions"""

def fTheta(x):
    return np.sin(x) / 2

def fphiPiPi(x):
    return (3 - np.cos(2 * x)) / (6 * np.pi)

def fphiJPsi(x):
    return (6 + np.cos(2 * x)) / (12 * np.pi)

"""2d functions"""

def f24(x2, x4):
    return (8 - 2*np.cos(2*x6) + np.cos(2*(-x2 + x4)) + np.cos(2 * (x2 + x4))) * np.sin(x2)

def f34(x3, x4):
    return (6 + np.cos(-2*x3) + 2 * np.cos(2*x4) + np.cos(2 * (x3 + x4)))

def f23(x2, x3):
    return (8 + 2*np.cos(-2*x3) * np.sin(x2) + cos(-2*x3)) * np.sin(3*x2)

def f14(x1, x4):
    return (8 - 2*np.cos(2*x4) + 2 * np.cos(2*x4) * np.cos(2*x1)) * np.sin(x1)

def f12(x1, x2):
    return ((3 + np.cos(2 * x1))*np.cos(x2)*np.cos(x2) - (3 + np.cos(x1)*np.cos(x1))*np.sin(x2)*np.sin(x2))

