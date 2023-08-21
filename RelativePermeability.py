# -*- coding: utf-8 -*-
"""
Created on Fri May 27 07:45:21 2022

@author: Maryelin
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

def relative_permeability(kroiw, krwro, krgro, swir, sorw, sorg, sgc,
                          nw, nsorw, ng, nog):
    #variables
    sw = np.linspace(swir, 1 - sorw, 20, endpoint=True)
    sg = np.linspace(sgc, 1 - sorg, 20, endpoint=True)
    so = 1 - sg
    #Models Corey, 1954
    krw = krwro * ((sw - swir) / (1 - sorw - swir))**nw

    krow = kroiw * ((1 - sw - sorw) / (1 - sorw - swir))**nsorw

    krg = krgro * ((sg - sgc) / (1 - sgc - sorg - swir))**ng
    krg[krg >=1] = 1

    krog = kroiw * ((1 - sg - sorg - swir) / (1 - sgc - sorg - swir))**nog

    #Stone's Models II

    kro = kroiw * ((((krow / kroiw) + krw) * ((krog / kroiw) + krg)) - krw - krg)
    kro[kro < 0] = 0 #in the case of negative values, KRO[i<0] = 0 Bakker, 2015
    return sw, krw, sg, krg, so, kro

# sw, krw, sg, krg, so, kro = relative_permeability(kroiw= 1, krwro=1, krgro=1, swir=0.1, sorw=0.1, sorg=0.1, sgc=0.1, nw=2.0, nsorw=2.0, ng=2.0, nog=2.0)



def stone_model_I(swir, sorg, sorw, sgc, krwro, kroiw, krgro, nw, nsorw, ng, nog):
    assert swir < 1
    #oil-water system and gas-oil system Corey two phases model
    #variables
    sw = np.linspace(swir, 1 - sorw, 20, endpoint=True)
    sg = np.linspace(sgc, 1 - sorg, 20, endpoint=True)
    so = 1 - sg
    #Models Corey, 1954
    krw = krwro * ((sw - swir) / (1 - sorw - swir))**nw

    krow = kroiw * ((1 - sw - sorw) / (1 - sorw - swir))**nsorw

    krg = krgro * ((sg - sgc) / (1 - sgc - sorg - swir))**ng
    krg[krg >=1] = 1

    krog = kroiw * ((1 - sg - sorg - swir) / (1 - sgc - sorg - swir))**nog

    #Stone Model I normalized by Aziz and Settari, 1979
    #swc = swir
    #Fayers and Mattews 1984
    a = 1 - (sg / (1 - swir - sorg))
    som= (a * sorw) + ((1 - a) * sorg)
    s_o = np.abs(so - som) / (1 - swir - som)  # so>= som
    s_w = np.abs(sw - swir) / (1 - swir - som)  # sw >= swir
    s_g = (sg) / (1 - swir - som)
    s_o[s_o >= 1.0] = 1 - swir
    s_w[s_w >= 1.0] = 1 - sorw
    s_g[s_g >= 1.0] = 1 - sorg
    kro0 = kroiw
    kro = (s_o / kro0) * (krow / (1 - s_w)) * (krog / (1 - s_g))
    kro[kro >= 1] = 1
    plt.plot(sw, kro,label='kro')
    plt.plot(sw, krw, label='krw')
    plt.ylabel('Relative Permability')
    plt.xlabel('sw')
    plt.grid(True, linewidth=1, which="both")
    plt.legend(loc='upper center')
    plt.show()
    return sw, krw, sg, krg, so, kro
sw, krw, sg, krg, so, kro = stone_model_I(swir=0.1, sorg=0.1, sorw=0.1, sgc=0.1, krwro=0.9, kroiw=1, krgro=0.9, nw=2, nsorw=2, ng=2, nog=2)

#
# plt.plot(sg, krg, label='krg')
# plt.plot(sg, krog, label='krog')
# plt.ylabel('Relative Permability')
# plt.xlabel('Sg')

# print(sw, krw, sg, krg, so, kro)
# s_w[s_w == 0] = s_w[np.abs(s_w - 0).argmin() +1]

def permeability(swir, sorg, sorw, sgc, krwro, kroiw, krgro, kroge, nw, nsorw ):
    assert swir < 1
    sw = np.linspace(swir, 1 - sorw, 20, endpoint=True)
    sg = np.linspace(sgc, 1 - sorg, 20, endpoint=True)
    so = 1 - sg

    #2D models Bakker 2015
    #Water-Oil=
    swd = (sw - swir) / (1 - swir - sorw)
    krw = krwro * (swd)**nw
    krow = kroiw * (1 - swd)**nsorw
    # print(krw)
    # print(krow)
    # plt.plot(sw, krw, 'r')
    # plt.plot(sw, krow, 'g')
    #Gas-Oil
    sgd = (sg) / (1 - swir - sorg)
    sgcd = (sg - sgc) / (1 - sgc - swir - sorg)
    krg = krgro * (sgd)**2 * (sgcd)**2
    krog = kroge * (1 - sgd)**4
    print(krg)
    print(krog)
    plt.plot(sg, krg)
    plt.plot(sg,krog)
# permeability(swir=0.1, sorg=0.1, sorw=0.1, sgc=0.1, krwro=1, kroiw=1, krgro=1, kroge=1, nw=2, nsorw=2)