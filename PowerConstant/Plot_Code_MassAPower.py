# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:10:07 2022

@author: Maryelin
"""
import itertools
import os
import numpy as np
import re
import pandas as pd
import itertools
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
import matplotlib.lines as Line2D

"Production"

# path3000 = 'production_Kerogen_Power_Cta.py_3000.txt'
# path2500 = 'production_Kerogen_Power_Cta.py_2500.txt'
# path2000 = 'production_Kerogen_Power_Cta.py_2000.txt'
# path1500 = 'production_Kerogen_Power_Cta.py_1500.txt'
# files = [path3000, path2500, path2000, path1500]

# color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in files:
#     pathfile = os.path.join(os.getcwd(), file)
#     with open(pathfile, 'r') as f:
#         data = pd.read_csv(pathfile, sep=' ', engine='python', header=None)
#         time = np.linspace(0, (data.iloc[:,1].max()),len(data))
#         c = next(color)
#         label = re.findall(r'\d+', file)[0]
#         Temp = data.iloc[:,6].max()
#         # ax.plot(time, data.iloc[:,2]/1000, c=c, label=f'{label} Watts Tmax={round(Temp,0)} K')
#         # yhat = savgol_filter(data.iloc[:,2]/100, 51, 3)
#         # ax.plot(time, yhat, c=c, label=f'{label}' ' Watts $\mathregular{T_{max}}$='f'{round(Temp, 0)} K')
#         ax.plot(time, data.iloc[:,5], c=c, label=f'{label} Watts')

# # ax.set_title('(a)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# # ax.set_ylabel('cumulative mass of oil produced ( $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Temperature (K)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper left', fontsize = '8' )

# plotpath= f'power_ctte_temperature.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')

# plt.show()


"KEROGEN"
'3000 W'
# pathfolder  = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Power_Cta.py_3000')
# files = os.listdir(pathfolder)
# r = re.compile(r"(\d+)")
# files.sort(key=lambda x: int(r.search(x).group(1)))

# x = []
# path = os.path.join(pathfolder, 'step0.txt')
# with open (path, 'r') as f:
#     read = f.readlines()[0:2500]
#     total_kero_ini = 0
#     for line in read:
#         lines = line.split()
#         total_kero_ini += (float(lines[3]))
#         x.append(float(lines[2]))
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
#             mass_kero += (float(content[3]))
#     total_kero.append(mass_kero)
# total_oil = []
# for file in files:
#     path = os.path.join(pathfolder, file)
#     with open(path, 'r') as f:
#         read = f.readlines()[0:2501]
#         massoil = 0
#         for line in read:
#             content = line.split()
#             massoil += float(content[4])
#         total_oil.append(massoil)
# '2500'
# pathfolder  = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Power_Cta.py_2500')
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
# '2000 W'
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Power_Cta.py_2000')
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

# '1500 W'
# pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Power_Cta.py_1500')
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

# 'Plot Power Kerogen'
# time = np.linspace(0, np.max(time), len(time))
# time1 = np.linspace(0, np.max(time1), len(time1))
# time2 = np.linspace(0, np.max(time2), len(time2))
# time3 = np.linspace(0, np.max(time3), len(time3))

# fig, ax = plt.subplots()
# ax = plt.subplot()

# ax.plot(time,[i/1000000 for i in total_kero],'tab:blue', label = '3000 Watts')
# ax.plot(time1,[i/1000000 for i in total_kero1],'tab:orange', label = '2500 Watts')
# ax.plot(time2,[i/1000000 for i in total_kero2],'tab:green', label = '2000 Watts')
# ax.plot(time3,[i/1000000 for i in total_kero3],'tab:red', label = '1500 Watts')
# # ax1 = ax.twinx()
# # ax1.plot(time, [i/1000000 for i in total_oil],'tab:blue', linestyle = '--')
# # ax1.plot(time1, [i/1000000 for i in total_oil1],'tab:orange', linestyle = '--')
# # ax1.plot(time2, [i/1000000 for i in total_oil2],'tab:green', linestyle = '--')
# # ax1.plot(time3, [i/1000000 for i in total_oil3],'tab:red', linestyle = '--')

# k = np.linspace(2.5077067493855306, 2.968312500000092, 6)
# k2 = np.linspace(30.279854107593004, 30.72977967320378, 6)

# ax.set_title('(b)', loc = 'left')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_yticks(k)
# ax1.set_ylabel('Total Mass of Oil (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax1.set_yticks(k2)
#
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper left', fontsize = '8', bbox_to_anchor=(0, 0.60))

# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Kerogen'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Oil')]
# plt.legend(handles=legend2, loc='center left', fontsize = '8', bbox_to_anchor=(0, 0.35))
# #
# plotpath= f'powewr_ctte_kerogen.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#
# plt.show()
"Boundary Comparison"
# path3000 = 'production_Kerogen_Power_Cta.py_3000.txt'
# path2500 = 'production_Kerogen_Power_Cta.py_2500.txt'
# path2000 = 'production_Kerogen_Power_Cta.py_2000.txt'
# path1500 = 'production_Kerogen_Power_Cta.py_1500.txt'
# # temp_1200 = 'production_Temp_Cta.py_1200.txt'
# # temp_1000 = 'production_Temp_Cta.py_1000.txt'
# # temp_900 = 'production_Temp_Cta.py_900.txt'
# # temp_800 = 'production_Temp_Cta.py_800.txt'
# files = [path3000, path2500, path2000, path1500]
# # files_ = [temp_1000, temp_800]

# color = ['tab:blue', 'tab:orange']
# color = itertools.cycle(color)
# fig, ax = plt.subplots()
# for file in files:
#     pathfile = os.path.join(os.getcwd(), file)
#     with open(pathfile, 'r') as f:
#         data = pd.read_csv(pathfile, sep=' ', engine='python', header=None)
#         time = np.linspace(0, (data.iloc[:,1].max()),len(data))
#         c = next(color)
#         label = re.findall(r'\d+', file)[0]
#         Temp = data.iloc[:, 6].max()
#         ax.plot(time, data.iloc[:, 2] / 1000, c=c, label=f'{label} Watts Tmax={round(Temp,0)} K')
#         yhat = savgol_filter(data.iloc[:, 2] / 1000, 51, 3)
#         ax.plot(time, yhat, c=c, label=f'{label}' ' Watts $\mathregular{T_{max}}$='f'{round(Temp, 0)} K')


# for file in files_:
#     pathfile = os.path.join(os.getcwd(), file)
#     with open(pathfile, 'r') as f:
#         data = pd.read_csv(pathfile, sep=' ', engine='python', header=None)
#         time = np.linspace(0, (data.iloc[:,1].max()),len(data))
#         c = next(color)
#         label = re.findall(r'\d+', file)[0]
#         # ax.plot(time, data.iloc[:,2]/1000, c=c, linestyle='--', label=f'{label} K')
#         yhat = savgol_filter(data.iloc[:, 2] / 1000, 51, 3)
#         ax.plot(time, yhat, c=c, linestyle='--', label=f'{label} K')
# ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_ylabel('cumulative mass of oil produced (x $\mathregular{10^{3}}$ Kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper left', fontsize = '8' )

# plotpath= f'comparison_boundaries.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')

# plt.show()
#  plt.close()

"Kerogen-both"

'2000 W'
pathfolder = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Power_Cta.py_2000')
files = os.listdir(pathfolder)
r = re.compile(r"(\d+)")
files.sort(key=lambda x: int(r.search(x).group(1)))

step2 = []
time2 = []
for file in files:
    path = os.path.join(pathfolder, file)
    with open(path, 'r') as f:
        readline = f.readline().split()
        step2.append(int(float(readline[0])))
        time2.append(float(readline[1]))
total_kero2 = []
for file in files:
    path = os.path.join(pathfolder, file)
    with open(path, 'r') as f:
        read = f.readlines()[0:2500]
        mass_kero = 0
        for line in read:
            content = line.split()
            mass_kero += (float(content[3]))
        total_kero2.append(mass_kero)
total_oil2 = []
for file in files:
    path = os.path.join(pathfolder, file)
    with open(path, 'r') as f:
        read = f.readlines()[0:2501]
        massoil = 0
        for line in read:
            content = line.split()
            massoil += float(content[4])
        total_oil2.append(massoil)
        
'Temperature'        
pathfolder  = os.path.join(os.getcwd(), 'Model_Results_Kerogen_Temp_Cta.py_800')
files = os.listdir(pathfolder)
r = re.compile(r"(\d+)")
files.sort(key=lambda x: int(r.search(x).group(1)))
step1 = []
time1 = []
for file in files:
    path = os.path.join(pathfolder, file)
    with open (path, 'r') as f:
        readline = f.readline().split()
        step1.append(int(float(readline[0])))
        time1.append(float(readline[1]))
total_kero1 = []
for file in files:
    path = os.path.join(pathfolder, file)
    with open (path, 'r') as f:
        read = f.readlines()[0:2500]
        mass_kero = 0
        for line in read:
            content = line.split()
            mass_kero += (float(content[3]))
        total_kero1.append(mass_kero)
        
time2000 = np.linspace(0, np.max(time2), len(time2))
time800 = np.linspace(0, np.max(time1), len(time1))

fig, ax = plt.subplots()
ax = plt.subplot()
        
ax.plot(time2000, [i/1000000 for i in total_kero2],'tab:green', label = '2000 Watts')
ax.plot(time800,[i/1000000 for i in total_kero1],'tab:red', label = '800 K')

ax.set_xlabel('Time (Year)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.set_ylabel('Total Mass of Kerogen (x $\mathregular{10^{6}}$ kg)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')


ax.grid(True, linewidth=1, which="both")
ax.legend(loc='upper left', fontsize = '8', bbox_to_anchor=(0, 0.60))

# #
# plotpath= f'powewr_ctte_kerogen.png'
# plt.savefig(plotpath, format='png', bbox_inches = 'tight')
#
plt.show()

