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
from matplotlib.pyplot import cm
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
import itertools
import collections

"BENCHMARK MODEL"
'Mass and Volumen Produccion, Viscosity'
prod_in_situ = 'BM_production.txt'
prod_NoKerogen = 'NoKerogen.py.txt'
prod = [prod_in_situ, prod_NoKerogen]
color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
color = itertools.cycle(color)
fig, ax = plt.subplots()
for file in prod:
    path = os.path.join(os.getcwd(), file)
    with open(path, 'r') as file:
        data = pd.read_csv(file, sep= ' ', engine='python', header=None)
        time = np.linspace(0, data.iloc[:,1].max(),len(data))
        c = next(color)
        oil = np.linspace(np.min((data.iloc[:,2] - data.iloc[0,2])/1000), np.max((data.iloc[:,2] - data.iloc[0,2])/1000), 10)
        # ax.plot(data.iloc[:,1], (data.iloc[:,2] - data.iloc[0,2])/1000, c=c)
        yhat = savgol_filter((data.iloc[:, 2]) * 0.00110231, 51, 3)
        ax.plot(time, yhat, c=c)
ax.set_title('(a)', loc = 'left')
ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced ($\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.set_ylabel('cumulative mass of oil produced (ton)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.grid(True, linewidth=1, which="both")
legend1 = [plt.Line2D([0],[0], color='tab:blue', label='In Situ Kerogen Decomposition'),
            plt.Line2D([0],[0], color='tab:orange', label='Primary Recovery')]
plt.legend(handles=legend1, loc='upper left', fontsize = '8')

# plotpath= f'figure7.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#


# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='Oil')]
# plt.legend(handles=legend1, loc='center right', fontsize = '8')
plt.show()

'Temperature in Cell near of injection'
# prod_in_situ = 'BM_production.txt'
# prod = [prod_in_situ]
# color = ['tab:blue']
# fig, ax = plt.subplots()
# with open(prod_in_situ, 'r') as file:
#     data = pd.read_csv(file,sep=' ', engine='python', header=None)
#     time = np.linspace(0, data.iloc[:, 1].max(), len(data))
#     ax.plot(time,data.iloc[:,12],'tab:orange', label='Rock')
#     ax.plot(time,data.iloc[:,11], 'tab:green', label='Kerogen')
#     ax.plot(time,data.iloc[:,13], 'black', linestyle='--',label='Temperature decomposition')
# # ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Temperature (K)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='lower right', fontsize = '8' )
# plt.show()
"profile temp"
# folderpath = os.path.join(os.getcwd(),'BM_Model_Results_temp_profile') # temperature
# files = os.listdir(folderpath)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in files:
#     path = os.path.join(folderpath, file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         step = data.iloc[0,0]
#         if step in zip([392, 1523]):
#             group = data.groupby(data.iloc[:2500,2]).mean()
#             group2 = data.groupby(data.iloc[:600,2]).mean() # [:601,2] to take the vistual cell
#             ax.plot(group[2], group[8]) #rock
#             ax.plot(group2[2], group2[7], linestyle='--') #kerogen
# # ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('x (m)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Temperature (K)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
#
# # legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Rock'),
# #           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Kerogen')]
# # plt.legend(handles=legend2, loc='center left', fontsize = '8')
# plt.show()
"profile pessure"
#
# folderpath1 = os.path.join(os.getcwd(),'Model_Results_pressure_profile(1mD)') # pressure
# files1 = os.listdir(folderpath1)
# r = re.compile(r"(\d+)")
# files1.sort(key=lambda x: int(r.search(x).group(1)))
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# # fig, ax = plt.subplots()
# for file1 in files1:
#     path1 = os.path.join(folderpath1, file1)
#     with open(path1, 'r') as file1:
#         data1 = pd.read_csv(file1, sep= ' ', engine='python', header=None)
#         step1 = data1.iloc[0,0]
#         if step1 in zip([1523]):
#             group1 = data1.groupby(data1.iloc[:2500,2]).mean()
#             # ax.plot(group1[2], group1[9])  # pressure rock
#
# folderpath2 = os.path.join(os.getcwd(),'Model_Results_pressure_profile(5mD)') # pressure
# files2 = os.listdir(folderpath2)
# r = re.compile(r"(\d+)")
# files2.sort(key=lambda x: int(r.search(x).group(1)))
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# # fig, ax = plt.subplots()
# for file2 in files2:
#     path2 = os.path.join(folderpath2, file2)
#     with open(path2, 'r') as file2:
#         data2 = pd.read_csv(file2, sep= ' ', engine='python', header=None)
#         step2 = data2.iloc[0,0]
#         if step2 in zip([7345]):
#             group2 = data2.groupby(data2.iloc[:2500,2]).mean()
#             # ax.plot(group[2], group[9])  # pressure rock
#
# folderpath3 = os.path.join(os.getcwd(),'Model_Results_pressure_profile(10mD)') # pressure
# files3 = os.listdir(folderpath3)
# r = re.compile(r"(\d+)")
# files3.sort(key=lambda x: int(r.search(x).group(1)))
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# # fig, ax = plt.subplots()
# for file3 in files3:
#     path3 = os.path.join(folderpath3, file3)
#     with open(path3, 'r') as file3:
#         data3 = pd.read_csv(file3, sep= ' ', engine='python', header=None)
#         step3 = data3.iloc[0,0]
#         if step3 in zip([7369]):
#             group3 = data3.groupby(data3.iloc[:2500,2]).mean()
#             # ax.plot(group3[2], group3[9])  # pressure rock
#
# fig, ax = plt.subplots()
# ax.plot(group1[2], group1[9] / 1000000, 'tab:blue', label='1mD')  # pressure rock
# ax.plot(group2[2], group2[9] / 1000000, 'tab:orange', label='5mD')  # pressure rock
# ax.plot(group3[2], group3[9] / 1000000, 'tab:green', label='10mD')  # pressure rock
#
# ax.set_title('(b)', loc = 'left')
# p = np.linspace(np.min(group3[9] / 1000000), np.max(group1[9] / 1000000), 5)
# ax.set_xlabel('x (m)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Pressure (MPa)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center right', fontsize = '8' )
# ax.set_yticks(p)
# plt.show()

'Kerogen'
# folderpath = os.path.join(os.getcwd(),'BM_Model_Results')
# files = os.listdir(folderpath)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# time = []
# kerogen = []
# oil = []
# for file in files:
#     path = os.path.join(folderpath, file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time.append(data.iloc[0,1])
#         kerogen.append((data.iloc[:2500,5].sum())/1000000)
#         # oil.append((data.iloc[:,6].sum())/1000000)
#
# ax.plot(time,kerogen,'tab:blue', label='Kerogen')
# # yhat = savgol_filter(kerogen, 51, 3)
# # ax.plot(time, yhat)
# k = np.linspace(np.min(kerogen), np.max(kerogen), 5)
# # ax1 = ax.twinx()
# # ax1.plot(time, oil, 'tab:blue', linestyle ='--', label='Oil')
# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen ($\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_yticks(k)
# # ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize='8')
# # legend2 = [plt.Line2D([0],[0], color = 'tab:blue', linestyle='-', label='Kerogen'),
# #           plt.Line2D([0],[0], color = 'tab:blue', linestyle='--', label='Oil')]
# # plt.legend(handles=legend2, loc='center left', fontsize = '8')
# plotpath= f'benchmark_model_kerogen_decomp.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()
'kerogen  decomposition'
# prod_in_situ = 'production_kerogen.py_900.txt'
# ax = plt.subplot()
# with open(prod_in_situ, 'r') as file:
#     data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#     time = np.linspace(0, data.iloc[:, 1].max(),len(data))
#     ax.plot(time, (data.iloc[0, 6] - data.iloc[:, 6])/1000, 'tab:blue', label='Oil Released by decomposition')
#     ax.plot(time, (data.iloc[:, 7] - data.iloc[0, 7])/1000, 'tab:orange', label='Oil Remaining In Rock Porous')
#     ax.plot(time, (data.iloc[:,2] - data.iloc[0,2])/1000,'tab:green', label='Oil Produced')
#     # ax.plot(time, ((data.iloc[:, 7] - data.iloc[0, 7])/1000) + ((data.iloc[:,2] - data.iloc[0,2])/1000),'tab:red', label='test')
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Mass of oil (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper left', fontsize = '8' )
## plotpath= f'BM_comparison_oil.png'
## plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()
# "HEATER POSITION"
# 'Mass Production, Temperature'
# distance_5m = 'Production_Position5m.txt'
# distance_10m = 'Production_Position10m.txt'
# distance_20m = 'Production_Position20m.txt'
# prod = [distance_5m, distance_10m, distance_20m]
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in prod:
#     path = os.path.join(os.getcwd(), file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time = np.linspace(0, data.iloc[:,1].max(),len(data))
#         c = next(color)
#         ax.plot(time, data.iloc[:,2]/1000, c=c)
#         # ax.plot(time, data.iloc[:,9], c=c)

# # ax.set_title('(a)', loc = 'left')
# # ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Temperature (K)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")

# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='5 m'),
#           plt.Line2D([0],[0], color='tab:orange', label='10 m'),
#             plt.Line2D([0],[0], color='tab:green', label='20 m')]
# plt.legend(handles=legend1, loc='upper left', fontsize = '8')
# # plt.legend(handles=legend1, loc='center right', fontsize = '8')

# plotpath= f'figure_14b.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()
# plt.close()
'Kerogen'
# # folderpath_5m = os.path.join(os.getcwd(),'Model_Results(5m)')
# # folderpath_10 = os.path.join(os.getcwd(),'Model_Results(10m)')
# # folderpath_20 = os.path.join(os.getcwd(),'Model_Results(20m)')
# # Folders = [folderpath_5m, folderpath_20]
# # fig, ax = plt.subplots()
# # color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# # color = itertools.cycle(color)
# # for folder in Folders:
# #     files = os.listdir(folder)
# #     r = re.compile(r"(\d+)")
# #     files.sort(key=lambda x: int(r.search(x).group(1)))
# #     time = []
# #     kerogen = []
# #     oil = []
# #     for file in files:
# #         path = os.path.join(folder, file)
# #         with open(path, 'r') as file:
# #             data = pd.read_csv(file, sep= ' ', engine='python', header=None)
# #             time.append(data.iloc[0,1])
# #             kerogen.append((data.iloc[:2500,3].sum())/1000000)
# #             oil.append((data.iloc[:,4].sum())/1000000)
# #     c = next(color)
# #     ax.plot(time,kerogen,c=c, label='Kerogen')
# #     ax1 = ax.twinx()
# #     ax1.plot(time, oil, c=c, linestyle ='--', label='Oil')
#
"THERMAL PROPERTIES"
'Rock Heat Capacity(ca_mc)'
# ####MASS PRODUCTION##
# ca_mc500 = 'production_ca_mc(500).txt'
# ca_mc1000 = 'production_ca_mc(1000).txt' #Benchmark
# ca_mc1500 = 'production_ca_mc(1500).txt'
# ca_mc2000 = 'production_ca_mc(2000).txt'
# prod = [ca_mc500, ca_mc1000, ca_mc1500, ca_mc2000]
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
#         # ax.plot(data.iloc[:,1], data.iloc[:,2]/1000, c=c, label=f'{label} J/(kg x K)') #mass
#         # ax.plot(data.iloc[:,1], data.iloc[:, 6] / 1000000, c=c, label=f'{label} J/(kg x K)')
#         yhat = savgol_filter(data.iloc[:, 2] / 1000, 51, 3)
#         ax.plot(time, yhat, c=c, label=f'{label} J/(kg x K)')
# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
#
# plotpath= f'rock_heat_capacity_production.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#
# plt.show()
# plt.close()

#Kerogen##

# ca_mc500 = os.path.join(os.getcwd(),'Model_Results_ca_mc(500)')
# ca_mc1000 = os.path.join(os.getcwd(),'Model_Results_ca_mc(1000)')
# ca_mc1500 = os.path.join(os.getcwd(),'Model_Results_ca_mc(1500)')
# ca_mc2000 = os.path.join(os.getcwd(),'Model_Results_ca_mc(2000)')
# prod = [ca_mc500, ca_mc1000, ca_mc1500, ca_mc2000]
# fig, ax = plt.subplots()
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# for folder in prod:
#     files = os.listdir(folder)
#     r = re.compile(r"(\d+)")
#     files.sort(key=lambda x: int(r.search(x).group(1)))
#     label = re.findall(r"\d+", folder)[0]
#     time = []
#     kerogen = []
#     oil = []
#     for file in files:
#         path = os.path.join(folder, file)
#         with open(path, 'r') as file:
#             data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#             time.append(data.iloc[0,1])
#             kerogen.append((data.iloc[:2500,5].sum())/1000000)
#             oil.append((data.iloc[:,6].sum())/1000000)
#     c = next(color)
#     ax.plot(time,kerogen,c=c, label=f'{label} J/(kg x K)')
    # ax1 = ax.twinx()
    # ax1.plot(time, oil, c=c, linestyle ='--', label='Oil')

# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# # ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# #
# # plotpath= f'thermal_properties_half_life_time_kerogen_1.png'
# # plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# #
# plt.show()
# # #  plt.close()

'Heat Transfer coefficient (ca_g)'
#MASS PRODUCTION##
# ca_g50 = 'production_ca_g(10).txt'
# ca_g100 = 'production_ca_g(100).txt'
# ca_g150 = 'production_ca_g(150).txt'
# prod = [ca_g50, ca_g100, ca_g150]
#
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
#         ax.plot(data.iloc[:,1], data.iloc[:,2]/1000, c=c, label=f'{label}' 'W/($\mathregular{m^{2}}$ x K)') #mass
        # ax.plot(time, data.iloc[:, 6] / 1000000, c=c, label=f'{label}' 'W/($\mathregular{m^{2}}$ x K)') #kerogen
# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )

# plotpath= f'rock_heat_capacity_production.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')

# plt.show()
# plt.close()

##Kerogen##
# folderpath_2 = os.path.join(os.getcwd(),'Model_Results_ca_g(50)')
# folderpath_3 = os.path.join(os.getcwd(),'Model_Results_ca_g(100)')
# folderpath_4 = os.path.join(os.getcwd(),'Model_Results_ca_g(150)')
# Folders = [folderpath_2, folderpath_3, folderpath_4]
#
# fig, ax = plt.subplots()
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
#
# for folder in Folders:
#     files = os.listdir(folder)
#     r = re.compile(r"(\d+)")
#     files.sort(key=lambda x: int(r.search(x).group(1)))
#
#     time = []
#     kerogen = []
#     oil = []
#     for file in files:
#         path = os.path.join(folder, file)
#         with open(path, 'r') as file:
#             data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#             time.append(data.iloc[0,1])
#             kerogen.append((data.iloc[:2500,5].sum())/1000000)
#             oil.append((data.iloc[:,6].sum())/1000000)
#     c = next(color)
#     label = re.findall(r"\d+", folder)[0]
#     ax.plot(time,kerogen,c=c, label=f'{label}' 'W/($\mathregular{m^{2}}$ x K)')
#     ax1 = ax.twinx()
#     ax1.plot(time, oil, c=c, linestyle ='--', label='Oil')
# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )

# plotpath= f'thermal_properties_half_life_time_kerogen_1.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')

# plt.show()
#  plt.close()

'Rock Heat Conductivity (fa_heatg)'
# ##MASS PRODUCTION ##
# # fa_heat1 = 'production_fa_heatg(1).txt'
# # fa_heat5 = 'production_fa_heatg(5).txt'
# # fa_heat15 = 'production_fa_heatg(15).txt'
# # prod = [fa_heat15, fa_heat10, fa_heat5, fa_heat1]
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
#         ax.plot(time, data.iloc[:,2]/1000, c=c, label=f'{label} W/(m x K)') #mass
#         # ax.plot(time, data.iloc[:, 6] / 1000000, c=c, label=f'{label} W/(m x K)')
# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# ax.legend(loc='center right', fontsize = '8' )
#
# # plotpath= f'rock_heat_capacity_production.png'
# # plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#
# plt.show()

'Thermal properties comparison'
# ca_mc500 = 'production_ca_mc(500).txt'
# ca_g2 = 'production_ca_g(150).txt'
# fa_heat15 = 'production_fa_heatg(15).txt'
# files = [fa_heat15, ca_mc500, ca_g2]
# # files = [fa_heat15, ca_mc500]
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in files:
#     path = os.path.join(os.getcwd(), file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time = np.linspace(0, data.iloc[:,1].max(),len(data))
#         c = next(color)
#         ax.plot(time, data.iloc[:,2]/1000, c=c) #mass
#         # ax.plot(time, data.iloc[:, 6] / 1000000, c=c)
#
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='Rock Heat Conductivity (15 W/(m x K))'),
#           plt.Line2D([0],[0], color='tab:orange', label='Rock Heat Capacity (500 J/(kg x K))'),
#             plt.Line2D([0],[0], color='tab:green', label='Rock Heat Transfer Coefficient (150 W/($\mathregular{m^{2}}$ x K))')]
# plt.legend(handles=legend1, loc='upper left', fontsize = '8')
# #
# # legend1 = [plt.Line2D([0],[0], color='tab:blue', label='Rock Heat Conductivity (15 W/(m x K))'),
# #           plt.Line2D([0],[0], color='tab:orange', label='Rock Heat Capacity (500 J/(kg x K))')]
# # plt.legend(handles=legend1, loc='upper left', fontsize = '8')

# plt.show()

'Activation Energy'
##MASS PRODUCTION##
# ka_de161 = 'production_ka_dE(161).txt'
# ka_de180 = 'production_ka_dE(180).txt'
# ka_de200 = 'production_ka_dE(200).txt'
# ka_de240 = 'production_ka_dE(240).txt'
#
# prod = [ka_de161, ka_de180, ka_de200, ka_de240]
#
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
#         ax.plot(time, data.iloc[:,2]/1000, c=c, label=f'{label} KJ/(kgMol)') #mass
#         # ax.plot(time, data.iloc[:, 5] / 1000000, c=c, label=f'{label}.{label2} W/(m x K)')
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8' )
# plt.show()
#Kerogen##

# pathfolder  = os.path.join(os.getcwd(), 'BM_Model_Results')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
# step0 = []
# time0 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         readline = f.readline().split()
#         step0.append(int(float(readline[0])))
#         time0.append(float(readline[1]))
# total_kero0 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[5]))
#         total_kero0.append(mass_kero)
# total_oil0 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil0.append(massoil)
#
# pathfolder  = os.path.join(os.getcwd(), 'Model_Results_ka_dE(180)')
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
#             mass_kero += (float(content[5]))
#         total_kero1.append(mass_kero)
# total_oil1 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil1.append(massoil)
#
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_ka_dE(200)')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
#
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
#             mass_kero += (float(content[5]))
#         total_kero2.append(mass_kero)
# total_oil2 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil2.append(massoil)
#
#
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_ka_dE(240)')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))
#
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
#             mass_kero += (float(content[5]))
#         total_kero3.append(mass_kero)
# total_oil3 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil3.append(massoil)
# time1 = np.linspace(0, np.max(time1), len(time1))
# time2 = np.linspace(0, np.max(time2), len(time2))
# time3 = np.linspace(0, np.max(time3), len(time3))
#
# fig, ax = plt.subplots()
# ax = plt.subplot()
#
# ax.plot(time0,[i/1000000 for i in total_kero0],'tab:blue', label = '161 KJ/(kgmol)')
# ax.plot(time1,[i/1000000 for i in total_kero1],'tab:orange', label = '180 KJ/(kgmol)')
# ax.plot(time2,[i/1000000 for i in total_kero2],'tab:green', label = '200 KJ/(kgmol)')
# ax.plot(time3,[i/1000000 for i in total_kero3],'tab:red', label = '240 KJ/(kgmol)')
# ax1 = ax.twinx()
# ax1.plot(time0, [i/1000000 for i in total_oil0],'tab:blue', linestyle = '--', label='161 KJ/(kgmol)')
# ax1.plot(time1, [i/1000000 for i in total_oil1],'tab:orange', linestyle = '--', label='180 KJ/(kgmol)')
# ax1.plot(time2, [i/1000000 for i in total_oil2],'tab:green', linestyle = '--', label='200 KJ/(kgmol)')
# ax1.plot(time3, [i/1000000 for i in total_oil3],'tab:red', linestyle = '--', label='240 KJ/(kgmol)')
#
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8',bbox_to_anchor=(0, 0.50))
#
# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.30))
#
# # plotpath= f'activation_energy_kerogen.png'
# # plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()

"Half Time"
##MASS PRODUCTION##
# ka_de180 = 'production_Kerogen.py_900_ka_de(180).txt'
# ka_de200 = 'production_Kerogen.py_900_ka_de(200).txt'
# ka_de240 = 'production_Kerogen.py_900_ka_de(240).txt'
# prod = [ka_de180, ka_de200, ka_de240]
#
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# labels = ['Half Year', 'One Year', 'Two Years']
# label = itertools.cycle(labels)
# fig, ax = plt.subplots()
# for file in prod:
#     path = os.path.join(os.getcwd(), file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time = np.linspace(0, data.iloc[:,1].max(),len(data))
#         c = next(color)
#         l = next(label)
#         ax.plot(time, data.iloc[:,2]/1000, c=c, label=l) #mass
        # ax.plot(time, data.iloc[:, 5] / 1000000, c=c)

##Kerogen##

# pathfolder  = os.path.join(os.getcwd(), 'Model_Results_ca_k2o(6m)')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))

# step = []
# time = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         readline = f.readline().split()
#         step.append(int(float(readline[0])))
#         time.append(float(readline[1]))
# total_kero = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open (path, 'r') as f:
#         read = f.readlines()[0:2500]
#         mass_kero = 0
#         for line in read:
#             content = line.split()
#             mass_kero += (float(content[5]))
#     total_kero.append(mass_kero)
# total_oil = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil.append(massoil)
# '1 year'
# pathfolder  = os.path.join(os.getcwd(), 'BM_Model_Results')
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
#             mass_kero += (float(content[5]))
#         total_kero1.append(mass_kero)
# total_oil1 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil1.append(massoil)
# '2 years'
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_ca_k2o(2y)')
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
#             mass_kero += (float(content[5]))
#         total_kero2.append(mass_kero)
# total_oil2 = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[6])
#         total_oil2.append(massoil)

# time = np.linspace(0, np.max(time), len(time))
# time1 = np.linspace(0, np.max(time1), len(time1))
# time2 = np.linspace(0, np.max(time2), len(time2))

# ax = plt.subplot()

# ax.plot(time,[i/1000000 for i in total_kero],'tab:blue', linestyle = '-' , label='Half year')
# ax.plot(time1,[i/1000000 for i in total_kero1],'tab:orange', linestyle = '-', label='One year')
# ax.plot(time2,[i/1000000 for i in total_kero2],'tab:green', linestyle = '-', label='Two years')
# ax1 = ax.twinx()
# ax1.plot(time, [i/1000000 for i in total_oil],'tab:blue', linestyle = '--')
# ax1.plot(time1, [i/1000000 for i in total_oil1],'tab:orange', linestyle = '--')
# ax1.plot(time2, [i/1000000 for i in total_oil2],'tab:green', linestyle = '--')

# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='center left', fontsize = '8',bbox_to_anchor=(0, 0.50))

# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.30))

# plotpath= f'time_half_life_kerogen.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()

"Permeability"
# BM_model = 'BM_production.txt'
# perm_10 = 'production_Kerogen.py_900_perm(5).txt'
# perm_20 = 'production_Kerogen.py_900_perm(10).txt'
# prod = [BM_model, perm_10, perm_20]
# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in prod:
#     path = os.path.join(os.getcwd(), file)
#     with open(path, 'r') as file:
#         data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#         time = np.linspace(0, data.iloc[:,1].max(),len(data))
#         c = next(color)
#         # ax.plot(data.iloc[:,1], (data.iloc[:,2] - data.iloc[0,2])/1000, c=c, linewidth='1.5')
#         yhat = savgol_filter(data.iloc[:, 2] / 1000, 51, 3)
#         ax.plot(time, yhat, c=c)
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# # ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='1 mD'),
#             plt.Line2D([0],[0], color='tab:orange', label='5 mD'),
#            plt.Line2D([0],[0], color='tab:green', label='10 mD')]
# plt.legend(handles=legend1, loc='center left', fontsize = '8')
#
# plt.show()

# prod_in_situ = 'production_perm(10mD).txt'
# ax = plt.subplot()
# with open(prod_in_situ, 'r') as file:
#     data = pd.read_csv(file, sep= ' ', engine='python', header=None)
#     time = np.linspace(0, data.iloc[:, 1].max(),len(data))
#     ax.plot(time, (data.iloc[0, 6] - data.iloc[:, 6])/1000, 'tab:blue', label='Oil Released by decomposition')
#     ax.plot(time, (data.iloc[:, 7] - data.iloc[0, 7])/1000, 'tab:orange', label='Oil Remaining In Rock Porous')
#     ax.plot(time,(data.iloc[:,2] - data.iloc[0,2])/1000,'tab:green', label='Oil Produced')
# ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Mass of oil (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper left', fontsize = '8' )
# plotpath= f'BM_comparison_oil.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
# plt.show()

"Custom Plots"
# ax.set_title('(a)', loc = 'left')
# ax.set_title('(b)', loc = 'left')
# # # # ax.set_title('(c)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_xlabel('x (m)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
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
# # ax.legend(loc='center right', fontsize = '8' )
# ax1.legend(loc='center left', fontsize = '8',bbox_to_anchor=(0, 0.35))

# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='Half Year'),
#           plt.Line2D([0],[0], color='tab:orange', label='One Year'),
#             plt.Line2D([0],[0], color='tab:green', label='Two Years')]
# plt.legend(handles=legend1, loc='center left', fontsize = '8')

# legend1 = [plt.Line2D([0],[0], color='tab:blue', label='In-situ Upgrading'),
#           plt.Line2D([0],[0], color='tab:orange', label='Heat Injection Only'),
#             plt.Line2D([0],[0], color='tab:green', label='Primary Recovery')]
# plt.legend(handles=legend1, loc='center right', fontsize = '8')

# legend1 = [plt.Line2D([0],[0], color='tab:orange', label='Heat Injection Only')]
# plt.legend(handles=legend1, loc='upper left', fontsize = '8')

# legend1 = plt.legend(loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.35))
# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.30))
# plt.gca().add_artist(legend1)
#
# plotpath= f'thermal_properties_half_life_time_kerogen_1.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#
# plt.show()
# #  plt.close()
