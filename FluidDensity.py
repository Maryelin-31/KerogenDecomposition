import numpy as np
from matplotlib import pyplot as plt

#Aalto-Keskinen Model Pure compound
#compute liquid molar volume
def liquid_density(vm, w, tcri, temp, pre, Pc):
    a = -1.52816
    b = 1.43907
    c = -0.81446
    d = 0.190454
    e = -0.296123
    f = 0.386914
    g = 0.0427258
    h = -0.0480645
    a0 = 48.85416
    a1 = -1154.2977
    a2 = 790.09727
    a3 = -212.14413
    a4 = 93.4904
    b0 = 0.0264002
    b1 = 0.42711522
    b2 = 0.5
    c1 = 9.2892236
    c2 = 2.5103968
    c3 = 0.5939722
    c4 = 0.0010895002
    D = 1.0001
    E = 0.80329503
    Tr = temp / tcri #Kelvin
    Pr = pre / Pc #MPA
    # Ps = Ps/Pc
    Ps = Pr
    #Saturate Liquid density Hankinson and Thomson
    vro = 1 + a * (np.abs(1 - Tr))**(1/3) + b * np.abs((1 - Tr))**(2/3) + c * (1 - Tr) + d * np.abs((1 - Tr))**(4/3)  # 0.25< Tr < 0.95
    vrd = (e + f*Tr + g*(Tr)**2 + h*(Tr)**3) / (Tr - 1.00001) #0.25 < Tr < 1.0
    vs = vm * (vro * (1 - w * vrd)) #m^3*mol^-1
    A = a0 + a1*Tr + a2*Tr**3 + a3*Tr**6 + a4/Tr
    B = b0 + ((b1) / (b2 + w))
    C = c1 * (1-Tr)**c2 + (1 - (1 - Tr)**c2)*np.exp(c3 + c4 * (Pr - Ps))
    v = 0.001 * (1 / (vs * ((A + (((C)**np.abs(D - Tr)**B) * np.abs(Pr -Ps)**E )) / (A + (C * np.abs(Pr - Ps )**E))))) #m^3/mol^-1
    return v
# temp = np.linspace(100, 700)
# pre = np.linspace(1.0e6, 12.0e6)
# v = liquid_density(vm=1267.5 * 1.0e-6, w = 0.751, tcri = 791.32 , temp = temp, pre=pre, Pc=1.25e6)







