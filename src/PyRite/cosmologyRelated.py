import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM

cosmo = LambdaCDM(H0 = 70, Om0= 0.3, Ode0= 0.7)

def findZ(t,cosmo = cosmo, dt = 1, zMax = 20, maxNstep = 1e5):
    """
    Give cosmic time in Myr and get redshift
    :param t: cosmic time in Myr
    :param cosmo: cosmology model, LCDM 70, 0.3, 0.7 default
    :param dt: accuracy of the redshift in Myr, default 1
    :param zMax: max redshift to test, default 20
    :param maxNstep: max number of iterrations, default 1e5
    :return: Log(SFR) in solar masses/yr
    """
    z0 = 0
    difference = 1e10
    nstep = 0
    ni = 10
    while abs(difference)>dt and maxNstep>nstep:
        currentGuess = 0
        diff = []
        xPos = []

        for i in range(ni):
            nstep += 1
            currentGuess = z0 + (zMax-z0) *i/(ni-1)
            xPos.append(currentGuess)
            difference = cosmo.age(currentGuess).value * 1e3 - t
            diff.append(difference)


            if i>0 and diff[-2] * difference < 0:
                currentGuess = zMax
                zMax = z0 + (zMax-z0) *i/(ni-1)
                z0 = z0 + (currentGuess-z0) *(i-1)/(ni-1)
                break

    return zMax

def MeridaMS(logM, z):
    """
    Give mass and redshift, get SFR
    :param M: Log(Mass) in solar masses
    :param z: Redshift
    :return: Log(SFR) in solar masses/yr
    """
    zPaper = np.array([0.0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375,
                  0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125,
                  0.875, 0.9375, 1.0, 1.0625, 1.125, 1.1875,
                  1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 6, 8.0,
                  8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
                  13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0,
                  17.5, 18.0, 18.5, 19.0, 19.5, 20.0])
    alphaPaper = np.array([-8.15625000000003, -8.13460937500003, -8.11218750000003,
                      -8.08898437500003, -8.06500000000003, -8.04023437500003,
                      -8.01468750000003, -7.988359375000028, -7.961250000000027,
                      -7.933359375000027, -7.904687500000026, -7.8752343750000255,
                      -7.8450000000000255, -7.813984375000024, -7.782187500000024,
                      -7.749609375000023, -7.716250000000022, -7.682109375000021,
                      -7.647187500000021, -7.611484375000019, -7.57, -7.27, -6.87,
                      -6.47, -6.33, -7.03, -7.75, -8.45, -8.578571428571419,
                      -8.688571428571416, -8.779999999999987, -8.852857142857129, -8.90714285714284,
                      -8.942857142857125, -8.96, -8.96,-8.96,
                      -8.96,-8.96,-8.96,-8.96,-8.96,-8.96,-8.96,-8.96,
                      -8.96,-8.96,-8.96,-8.96,-8.96,-8.96,-8.96,])
    betaPaper = np.array([0.8636250000000004, 0.8627109375000006, 0.8617187500000006,
                     0.8606484375000009, 0.8595000000000009, 0.8582734375000011,
                     0.8569687500000012, 0.8555859375000012, 0.8541250000000014,
                     0.8525859375000013, 0.8509687500000014, 0.8492734375000015,
                     0.8475000000000015, 0.8456484375000015, 0.8437187500000016,
                     0.8417109375000016, 0.8396250000000016, 0.8374609375000016,
                     0.8352187500000016, 0.8328984375000016, 0.83, 0.81, 0.78,
                     0.75, 0.75, 0.85, 0.93, 1.02,
                     1.0395238095238082, 1.0578571428571413, 1.0749999999999982,
                     1.090952380952379, 1.1057142857142832, 1.1192857142857116,
                     1.1316666666666637, 1.1428571428571395, 1.152857142857139,
                     1.1616666666666626, 1.16928571428571, 1.1757142857142808,
                     1.1809523809523756, 1.184999999999994, 1.1878571428571365,
                     1.1895238095238025, 1.1899999999999924, 1.189285714285706,
                     1.1873809523809435, 1.1842857142857048, 1.17999999999999,
                     1.1745238095237989, 1.1678571428571312, 1.1599999999999877])

    alphaInterp = np.interp(z,
                            zPaper,
                            alphaPaper)
    betaInterp = np.interp(z,
                            zPaper,
                            betaPaper)
    logSFR = betaInterp*logM + alphaInterp
    return logSFR

def PopessoMS(logM, z, IMF = 'C'):
    """
    Give mass and redshift, get SFR
    :param logM: Mass in solar masses
    :param z: Redshift
    :param IMF: C for chabrier, S for Salpeter, K for Kroupa
    :return: Log(SFR) in solar masses/yr
    """
    a0 = 2.71
    a1 = -0.186
    a2 = 10.86
    a3 = -0.0729


    age = cosmo.age(z).value

    logSFR = a0 + a1*age -np.log10(1 + 1/(10**logM/10**(a2+a3*age)))
    if IMF == 'K':
        return logSFR
    logSFR = np.log10(10**logSFR/0.67)
    if IMF == 'S':
        return logSFR
    logSFR = np.log10(10 ** logSFR*0.63)
    if IMF == 'C':
        return logSFR

    raise ValueError('Wrong IMF type!')

def KoprowskiMS(logM, z, IMF = 'C'):
    """
    Give mass and redshift, get SFR
    :param M: Log(Mass) in solar masses
    :param z: Redshift
    :param IMF: C for chabrier, S for Salpeter, K for Kroupa
    :return: Log(SFR) in solar masses/yr
    """
    a1 = 2.97
    a2 = -2.62
    a3 = 0.59
    b1 = 11.38
    b2 = -1.14
    b3 = 0.5

    sfrMax = 10**(a1 + a2*np.exp(-a3*z))
    M0 = 10**(b1 + b2*np.exp(-b3*z))

    SFR = sfrMax/(1+M0/10**logM)
    logSFR = np.log10(SFR)
    if IMF == 'C':
        return logSFR
    ValueError('Wrong IMF type! Write it you lazy dong!')

def PacificiMS(logM, z, QGLim = False, IMF = 'C'):
    """
    Give mass and redshift, get SFR
    :param M: Log(Mass) in solar masses
    :param z: Redshift
    :param QGLim: if True, returns limit for quiescence
    :param IMF: C for chabrier, S for Salpeter, K for Kroupa
    :return: Log(SFR) in solar masses/yr
    """
    cosmicAge = cosmo.age(z).value
    if QGLim:
        tau = np.log10(0.2 / (cosmicAge * 1e9))
    else:
        tau = np.log10(1 / (cosmicAge * 1e9))
    logSFR = tau+logM

    if IMF == 'C':
        return logSFR
    ValueError('Wrong IMF type! Write it you lazy dong!')
