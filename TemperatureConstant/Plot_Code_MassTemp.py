# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:10:07 2022

@author: Maryelin
"""

import os
import numpy as np
import re
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.lines as Line2D
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
import itertools
import decimal
    
"PRODUCTION"

# temp_1200 ='production_Kerogen_Temp_Cta.py_1200.txt'
# temp_1000 ='production_Kerogen_Temp_Cta.py_1000.txt'
# temp_900 = 'production_Kerogen_Temp_Cta.py_900.txt'
# temp_800 = 'production_Kerogen_Temp_Cta.py_800.txt'
# prod = [temp_1200, temp_1000, temp_900, temp_800]

# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in prod:
#     path = os.path.join(os.getcwd(), file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time = np.linspace(0, data.iloc[:,1].max(),len(data))
#         c = next(color)
#         label = re.findall(r"\d+", path)[0]
#         # ax.plot(data.iloc[:,1], data.iloc[:,2]/1000, c=c, label=f'{label} K') #mass
#         yhat = savgol_filter(data.iloc[:, 2] / 1000, 51, 3)
#         ax.plot(time, yhat, c=c, label=f'{label} K')
#         # ax.plot(time, data.iloc[:, 5] / 1000000, c=c, label=f'{label} K')
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# # plotpath= f'temp_ctte_production.png'
# # plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()

"KEROGEN"

# pathfolder  = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Temp_Cta.py_800')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
# step1 = []
# time1 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         readline = f.readline().split()
#         step1.append(int(float(readline[0])))
#         time1.append(float(readline[1]))
# total_kero1 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[3]))
#         total_kero1.append(mass_kero)

# total_oil1 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[4])
#         total_oil1.append(massoil)

# pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Temp_Cta.py_900')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))

# step2 = []
# time2 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         readline = f.readline().split()
#         step2.append(int(float(readline[0])))
#         time2.append(float(readline[1]))
# total_kero2 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[3]))
#         total_kero2.append(mass_kero)
# total_oil2 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[4])
#         total_oil2.append(massoil)


# pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Temp_Cta.py_1000')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))

# step3 = []
# time3 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         readline = f.readline().split()
#         step3.append(int(float(readline[0])))
#         time3.append(float(readline[1]))
# total_kero3 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[3]))
#         total_kero3.append(mass_kero)
# total_oil3 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[4])
#         total_oil3.append(massoil)
# #
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Temp_Cta.py_1200')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))

# step4 = []
# time4 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         readline = f.readline().split()
#         step4.append(int(float(readline[0])))
#         time4.append(float(readline[1]))
# total_kero4 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[3]))
#         total_kero4.append(mass_kero)
# total_oil4 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[4])
#         total_oil4.append(massoil)

# time1 = np.linspace(0, np.max(time1), len(time1))
# time2 = np.linspace(0, np.max(time2), len(time2))
# time3 = np.linspace(0, np.max(time3), len(time3))
# time4 = np.linspace(0, np.max(time4), len(time4))

# fig, ax = plt.subplots()
# ax = plt.subplot()

# ax.plot(time4,[i/1000000 for i in total_kero4],'tab:blue', label = '1200 K')
# ax.plot(time3,[i/1000000 for i in total_kero3],'tab:orange', label = '1000 K')
# ax.plot(time2,[i/1000000 for i in total_kero2],'tab:green', label = '900 K')
# ax.plot(time1,[i/1000000 for i in total_kero1],'tab:red', label = '800 K')
# # ax1 = ax.twinx()
# # ax1.plot(time4, [i/1000000 for i in total_oil4],'tab:blue', linestyle = '--', label='1200 K')
# # ax1.plot(time3, [i/1000000 for i in total_oil3],'tab:orange', linestyle = '--', label='1000 K')
# # ax1.plot(time2, [i/1000000 for i in total_oil2],'tab:green', linestyle = '--', label='900 K')
# # ax1.plot(time1, [i/1000000 for i in total_oil1],'tab:red', linestyle = '--', label='800 K')

# # ax.set_title('(b)', loc = 'left')
# k = np.linspace(2.5077067493855306, 2.968312500000092, 6)
# k2 = np.linspace(30.279854107, 30.740459853, 6)

# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_yticks(k)
# # ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# # ax1.set_yticks(k2)
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8',bbox_to_anchor=(0, 0.50))

# # legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
# #           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# # legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen')]
# # plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.30))

# plotpath= f'temp_ctte_kerogen.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()

"ENERGY"
# path800 = 'energy_Kerogen_Temp_Cta.py_800(2).txt'
# with open(path800, 'r') as file:
#     data800 = pd.read_csv(path800, sep=' ', engine='python', header=None)
#     data800 = data800.groupby(data800.iloc[:,1]).mean()
#     data800 = data800[data800.iloc[:] != 0].dropna()
#     group800 = data800.groupby(data800.iloc[:,1]).mean()
#     energy1_800 = data800.iloc[0, 3] - data800.iloc[:, 3]
#     energy2_800 = data800.iloc[0, 5] - data800.iloc[:, 5]
#     total_energy800 = (energy1_800 + energy2_800) / 1.0e9

# path900 = 'energy_Kerogen_Temp_Cta.py_900(2).txt'
# with open(path900, 'r') as file:
#     data900 = pd.read_csv(path900, sep=' ', engine='python', header=None)
#     data900 = data900.groupby(data900.iloc[:,1]).mean()
#     data900 = data900[data900.iloc[:] != 0].dropna()
#     group900 = data900.groupby(data900.iloc[:,1]).mean()
#     energy1_900 = data900.iloc[0, 3] - data900.iloc[:, 3]
#     energy2_900 = data900.iloc[0, 5] - data900.iloc[:, 5]
#     total_energy900 = (energy1_900 + energy2_900) / 1.0e9

# path1000 = 'energy_Kerogen_Temp_Cta.py_1000(2).txt'
# with open(path1000, 'r') as file:
#     data1000 = pd.read_csv(path1000, sep=' ', engine='python', header=None)
#     data1000 = data1000.groupby(data1000.iloc[:,1]).mean()
#     data1000 = data1000[data1000.iloc[:] != 0].dropna()
#     group1000 = data1000.groupby(data1000.iloc[:,1]).mean()
#     energy1_1000 = data1000.iloc[0, 3] - data1000.iloc[:, 3]
#     energy2_1000 = data1000.iloc[0, 5] - data1000.iloc[:, 5]
#     total_energy1000 = (energy1_1000 + energy2_1000) / 1.0e9

# path1200 = 'energy_Kerogen_Temp_Cta.py_1200(2).txt'
# with open(path1200, 'r') as file:
#     data1200 = pd.read_csv(path1200, sep=' ', engine='python', header=None)
#     data1200 = data1200.groupby(data1200.iloc[:,1]).mean()
#     data1200 = data1200[data1200.iloc[:] != 0].dropna()
#     group1200 = data1200.groupby(data1200.iloc[:,1]).mean()
#     energy1_1200 = data1200.iloc[0, 3] - data1200.iloc[:, 3]
#     energy2_1200 = data1200.iloc[0, 5] - data1200.iloc[:, 5]
#     total_energy1200 = (energy1_1200 + energy2_1200) / 1.0e9

# fig, ax = plt.subplots()
# ax.plot(data1200.iloc[:,1], total_energy1200, 'tab:blue', label='1200 K')
# ax.plot(data1000.iloc[:,1], total_energy1000, 'tab:orange', label='1000 K')
# ax.plot(data900.iloc[:,1], total_energy900, 'tab:green', label='900 K')
# ax.plot(data800.iloc[:,1], total_energy800, 'tab:red', label='800 K')

# # SMOOTH DATA OF ENERGY
# yhat1200 = savgol_filter(total_energy1200,51,3)
# yhat1000 = savgol_filter(total_energy1000,51,3)
# yhat900 = savgol_filter(total_energy900,51,3)
# yhat800 = savgol_filter(total_energy800,51,3)

# ax.plot(data1200.iloc[:,1], yhat1200, 'tab:blue', label='1200 K')
# ax.plot(data1000.iloc[:,1], yhat1000, 'tab:orange', label='1000 K')
# ax.plot(data900.iloc[:,1], yhat900, 'tab:green', label='900 K')
# ax.plot(data800.iloc[:,1], yhat800, 'tab:red', label='800 K')
# ax.grid(True, linewidth=1, which="both")
"POWER"
path800 = 'energy_Kerogen_Temp_Cta.py_800(2).txt'
with open(path800, 'r') as file:
    data800 = pd.read_csv(path800, sep=' ', engine='python', header=None)
    data800 = data800.groupby(data800.iloc[:,1]).mean()
    data800 = data800[data800.iloc[:] != 0].dropna()
#
path900 = 'energy_Kerogen_Temp_Cta.py_900(2).txt'
with open(path900, 'r') as file:
    data900 = pd.read_csv(path900, sep=' ', engine='python', header=None)
    data900 = data900.groupby(data900.iloc[:,1]).mean()
    data900 = data900[data900.iloc[:] != 0].dropna()

path1000 = 'energy_Kerogen_Temp_Cta.py_1000(2).txt'
with open(path1000, 'r') as file:
    data1000 = pd.read_csv(path1000, sep=' ', engine='python', header=None)
    data1000 = data1000.groupby(data1000.iloc[:,1]).mean()
    data1000 = data1000[data1000.iloc[:] != 0].dropna()

path1200 = 'energy_Kerogen_Temp_Cta.py_1200(2).txt'
with open(path1200, 'r') as file:
    data1200 = pd.read_csv(path1200, sep=' ', engine='python', header=None)
    data1200 = data1200.groupby(data1200.iloc[:,1]).mean()
    data1200 = data1200[data1200.iloc[:] != 0].dropna()

fig, ax = plt.subplots()
# ax.plot(data1200.iloc[:,1], data1200.iloc[:,2], 'tab:blue', label='1200 K')
# ax.plot(data1000.iloc[:,1], data1000.iloc[:,2], 'tab:orange', label='1000 K')
# ax.plot(data900.iloc[:,1], data900.iloc[:,2], 'tab:green', label='900 K')
# ax.plot(data800.iloc[:,1], data800.iloc[:,2], 'tab:red', label='800 K')

# smooth data to Power
yhat800 = savgol_filter(data800.iloc[:,2],51,3)
yhat900 = savgol_filter(data900.iloc[:,2],51,3)
yhat1000 = savgol_filter(data1000.iloc[:,2],51,3)
yhat1200 = savgol_filter(data1200.iloc[:,2],51,3)
ax.plot(data1200.iloc[:,1], yhat1200, 'tab:blue', label='1200 K')
ax.plot(data1000.iloc[:,1], yhat1000, 'tab:orange', label='1000 K')
ax.plot(data900.iloc[:,1], yhat900, 'tab:green', label='900 K')
ax.plot(data800.iloc[:,1], yhat800, 'tab:red', label='800 K')

"Custom Plots"
# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# # # # ax.set_title('(c)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative volume of oil produced ($\mathregular{m^{3}}$)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Mass of oil (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Viscosity (Pa*s)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Temperature (K)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# ax.legend(loc='center right', fontsize = '8' )
# ax.legend(loc='lower right', fontsize = '8' )
# ax.legend(loc='lower left', fontsize = '8' )
# ax.legend(loc='upper left', fontsize = '8' )
# # # ax.legend(loc='center right', fontsize = '8' )
# ax1.legend(loc='center left', fontsize = '8',bbox_to_anchor=(0, 0.35))

# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='Half Year'),
#           plt.Line2D([0],[0], color='tab:orange', label='One Year'),
#            plt.Line2D([0],[0], color='tab:green', label='Two Years')]
# plt.legend(handles=legend1, loc='center right', fontsize = '8')

# legend1 = plt.legend(loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.35))
# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.35))
# plt.gca().add_artist(legend1)
#
# plotpath= f'boundary_condition_tempconst_kerogen.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')

# plt.show()
# #  plt.close()
'Custom Energy Plots'

# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# ax.set_title('(c)', loc = 'left')
ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.set_ylabel('Power Injection (J/s)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# # ax.set_ylabel('Energy Input (x$\mathregular{10^{9}}$ J)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.grid(True, linewidth=1, which="both")
ax.legend(loc='upper right', fontsize = '8' )
# # ax.legend(loc='lower right', fontsize = '8' )
# plt.show()
# plt.close()

plotpath= f'temp_ctte_energy.png'
plt.savefig(plotpath, format='png', bbox_inches = 'tight')
plt.show()
