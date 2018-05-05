##Rate Coefficients

import numpy as np

def k1(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-32.71396786375
          + 13.53655609057*log_T_eV
          - 5.739328757388*log_T_eV**2 
          + 1.563154982022*log_T_eV**3
          - 0.2877056004391*log_T_eV**4
          + 0.03482559773736999*log_T_eV**5
          - 0.00263197617559*log_T_eV**6
          + 0.0001119543953861*log_T_eV**7
          - 2.039149852002e-6*log_T_eV**8)
    return rv

def k2(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv  = np.exp(-28.61303380689232
            - 0.7241125657826851*log_T_eV
            - 0.02026044731984691*log_T_eV**2
            - 0.002380861877349834*log_T_eV**3
            - 0.0003212605213188796*log_T_eV**4
            - 0.00001421502914054107*log_T_eV**5
            + 4.989108920299513e-6*log_T_eV**6
            + 5.755614137575758e-7*log_T_eV**7
            - 1.856767039775261e-8*log_T_eV**8
            - 3.071135243196595e-9*log_T_eV**9)
    return rv

def k3(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-44.09864886561001
             + 23.91596563469*log_T_eV
             - 10.75323019821*log_T_eV**2
             + 3.058038757198*log_T_eV**3
             - 0.5685118909884001*log_T_eV**4
             + 0.06795391233790001*log_T_eV**5
             - 0.005009056101857001*log_T_eV**6
             + 0.0002067236157507*log_T_eV**7
             - 3.649161410833e-6*log_T_eV**8)
    return rv

def k4(T):
    T_eV = T / 11605.

    rv = 3.92e-13/T_eV**0.6353
    return rv    

def k5(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-68.71040990212001
             + 43.93347632635*log_T_eV
             - 18.48066993568*log_T_eV**2
             + 4.701626486759002*log_T_eV**3
             - 0.7692466334492*log_T_eV**4
             + 0.08113042097303*log_T_eV**5
             - 0.005324020628287001*log_T_eV**6
             + 0.0001975705312221*log_T_eV**7
             - 3.165581065665e-6*log_T_eV**8)
    return rv  

def k6(T):
    rv = 3.36e-10/np.sqrt(T)/(T/1.**3)**0.2/ (1.+(T/1.**6)**0.7)
    return rv 

def k7(T):
    rv = 3.0e-16 * (T/3.**2)**0.95 * np.exp(-T/9.32e3)
    return rv 

def k8(T):
    rv = 1.35e-9*(T**9.8493e-2 + 3.2852e-1
        * T**5.5610e-1 + 2.771e-7 * T**2.1826)/ (1. + 6.191e-3 * T**1.0461
        + 8.9712e-11 * T**3.0424
        + 3.2576e-14 * T**3.7741)
    return rv 

def k9(T):
    rv = 2.10e-20 * (T/30.0)**(-0.15)
    return rv 

def k10(T):
    rv = 6.0e-10
    return rv 

def k11(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-24.24914687
                + 3.40082444*log_T_eV
                - 3.89800396*log_T_eV**2
                + 2.04558782*log_T_eV**3
                - 0.541618285*log_T_eV**4
                + 8.41077503*(10**-2)*log_T_eV**5
                - 7.87902615*(10**-3)*log_T_eV**6
                + 4.13839842*(10**-4)*log_T_eV**7
                - 9.36345888*(10**-6)*log_T_eV**8)
    return rv

def k12(T):
    rv = 5.6*(10**-11)*(T**0.5)*np.exp(-102124/T)
    return rv

def k13(T):
    T_eV = T / 11605.
    rv = 1.067*(10**-10)*(T_eV**2.012)*np.exp((-4.463/T_eV)*((1+0.2472*T_eV)**3.512))
    return rv

def k14(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-18.01849334
                + 2.3608522*log_T_eV
                - 0.28274430*log_T_eV**2
                + 1.62331664*log_T_eV**3
                - 3.36501203*log_T_eV**4
                + 1.17832978*(10**-2)*log_T_eV**5
                - 1.65619470*(10**-3)*log_T_eV**6
                + 1.06827520*(10**-4)*log_T_eV**7
                - 2.63128581*(10**-6)*log_T_eV**8)
    return rv

def k15(T):
    T_eV = T / 11605.
    log_T_eV = np.log(T_eV)
    rv = np.exp(-20.37260896
                    + 1.13944933*log_T_eV
                    - 0.14210135*log_T_eV**2
                    + 8.4644554*(10**-3)*log_T_eV**3
                    - 1.4327641*(10**-3)*log_T_eV**4
                    + 2.0122503*(10**-4)*log_T_eV**5
                    + 8.6639632*(10**-5)*log_T_eV**6
                    - 2.5850097*(10**-5)*log_T_eV**7
                    + 2.4555012*(10**-6)*log_T_eV**8
                    - 8.0683825*(10**-8)*log_T_eV**9)

    return rv

def k16(T):
    rv = 7*(10**-8)*((T/100)**-0.5)
    return rv

def k17(T):
    T_eV = T / 11605.
    if T_eV < 1.719:
        rv = 1.e-8*T**(-0.4)
    else:
        rv = 4.0e-4*T**(-1.4)*np.exp(-15100./T)
    return rv  

def k18(T):
    if T < 617:    
        rv = 1.e-8
    else:
        rv = 1.32e-6 * T**(-0.76)
    return rv  

def k19(T):
    rv = 5.e-7*np.sqrt(100./T)
    return rv   
