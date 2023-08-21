# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:10:07 2022

@author: Maryelin
"""

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import itertools
import re


fig, ax = plt.subplots()
folderpath = os.path.join(os.getcwd(),'Final_Result.txt')
with open(folderpath, 'r') as file:
    data = pd.read_csv(file, sep=' ', engine='python', header=None)
    # ax.plot(data.iloc[:, 0], data.iloc[:, 1], 'o', ms = 3.0, label='Simulation Results' ) #injection
    ax.plot(data.iloc[:, 0], data.iloc[:, 2], 'o', ms = 1)
    ax.plot(data.iloc[:, 0], data.iloc[:, 3],'o', ms = 1)
    ax.plot(data.iloc[:, 0], data.iloc[:, 4],'o', ms = 1)
    ax.plot(data.iloc[:, 0], data.iloc[:, 5],'o', ms = 1)
    # position = [1.6, 3.2, 4.8, 6.4]
    # ax.plot(position, data.iloc[2756,2:], 'tab:blue', label = f't = {round(data.iloc[2756,0], 1)} min')
    # ax.plot(position, data.iloc[2880,2:], 'green', label = f't = {round(data.iloc[3222,0], 1)} min')
    # ax.plot(position, data.iloc[2996,2:], 'red', label = f't = {round(data.iloc[3498,0], 1)} min')

# fig, ax = plt.subplots()
xscpath = os.path.join(os.getcwd(),'experiment.xlsx')
excel = pd.read_excel('experiment.xlsx', sheet_name='data')
# ax.plot(excel.iloc[:,1]  ,excel.iloc[:,2], '--',  linewidth=2, label='Experiment Results') #injection
ax.plot(excel.iloc[:,1] ,excel.iloc[:,3], 'tab:blue',linestyle='--', linewidth=1)
ax.plot(excel.iloc[:,1] ,excel.iloc[:,4], 'tab:orange', linestyle='--', linewidth=1)
ax.plot(excel.iloc[:,1] ,excel.iloc[:,5], 'tab:green', linestyle='--', linewidth=1)
ax.plot(excel.iloc[:,1] ,excel.iloc[:,6], 'tab:red', linestyle='--', linewidth=1)
# position = [1.6, 3.2, 4.8, 6.4]
# ax.plot(position ,excel.iloc[102,3:7], 'tab:blue', linestyle='--', linewidth=1, label = f't = {excel.iloc[102,1]} min')
# ax.plot(position ,excel.iloc[112,3:7], 'tab:orange', linestyle='--', linewidth=1)
# ax.plot(position ,excel.iloc[120,3:7], 'tab:green', linestyle='--', linewidth=1)
# temp = []
# t = []
# for row, col in excel.iterrows():
#     temp.append(col.iloc[3])
#     temp.append(col.iloc[4])
#     temp.append(col.iloc[5])
#     temp.append(col.iloc[6])
#     t.append(col.iloc[1])
# position = [1.6, 3.2, 4.8, 6.4]
# temp = np.reshape(temp, (377,4))
# temp = np.transpose(temp)
# plot = plt.contourf(temp,100, cmap='coolwarm', extent=[np.min(t), np.max(t), np.min(position), np.max(position)])


# fig, ax = plt.subplots()
# folderpath = os.path.join(os.getcwd(),'production_Kerogen_Power_Cta11boundary.py_1000.0.txt')
# temperature = []
# time = []
# with open(folderpath, 'r') as file:
#     readlines = file.readlines()
#     for line in readlines:
#         content = line.split()
#         temperature.append(float(content[2]))
#         temperature.append(float(content[3]))
#         temperature.append(float(content[4]))
#         temperature.append(float(content[5]))
#         time.append(int(float(content[0])))

# position = [1.6, 3.2, 4.8, 6.4]
# temperature = np.reshape(temperature, (8520,4))
# temperature = np.transpose(temperature)
# plot = plt.contourf(temperature,100, cmap='coolwarm', extent=[np.min(time), np.max(time), np.min(position), np.max(position)])

# plt.yticks(position)
# ax.set_xticks(position)
ax.set_xlabel('Time (min)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.set_xlabel('Monitoring Points Location (cm)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# plt.ylabel('Monitoring Points Location (cm)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
ax.set_ylabel('Temperature (C)', fontstyle='italic', fontname='Times New Roman', fontsize='12', fontweight='bold')
# ax.grid(True, linewidth=1, which="both")
# ax.legend(loc='upper right', fontsize = '8' )
# legend1 = plt.legend(loc='upper center', fontsize = '8')
# legend2 = [plt.Line2D([0],[0], color = 'k', linestyle='-', label='Simulation Results'),
#           plt.Line2D([0],[0], color = 'k', linestyle='--', label='Experiment Results')]
# plt.legend(handles=legend2, loc='upper right', fontsize = '8')
# plt.gca().add_artist(legend1)
# plt.clim(vmin=np.min(temperature), vmax=np.max(temperature))
# plt.clim(vmin=np.min(temp), vmax=np.max(temp))
# line = np.linspace(24.2, 31.2, 9)
# plt.clim(24.2, 31.2)
# plt.colorbar(ticks=line)
# plt.colorbar()
plt.show()



