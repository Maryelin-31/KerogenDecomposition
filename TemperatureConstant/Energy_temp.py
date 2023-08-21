import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as pt
from scipy import interpolate

# time=[]
# Power1=[]
# Power2=[]
# Energy1=[]
# Energy2=[]
# Temp_heater1=[]
# Temp_heater2=[]
# with open('energy.txt', 'r') as file:
#     for line in file:
#         read = line.split()
#         time.append(float(read[0]))
#         Power1.append(float(read[1]))
#         Power2.append(float(read[2]))
#         Energy1.append(float(read[3]))
#         Energy2.append(float(read[4]))
#         Temp_heater1.append(float(read[5]))
#         Temp_heater2.append(float(read[6]))
#
# Power = np.add(Power1, Power2) #Heat rate input
# Energy = np.add(Energy1, Energy2) #Energy Input
# print(Temp_heater1)
# #Thermal Accumulation Energy
# energy_1 = [] #Heater 1
# energini = Energy1[0]
# for i in Energy1:
#     energy_1.append(energini - i)
# energy_2 = [] #Heater 2
# energini = Energy2[0]
# for i in Energy2:
#     energy_2.append(energini - i)
# total_energy_accu = np.add(energy_1, energy_2)
#
# fig, ax1 = plt.subplots()
# ax1.plot(time, total_energy_accu, 'green', aa=True)

# ax2 = ax1.twinx()
# ax2.plot(time, Temp_heater1, 'red', aa=True)
#
# ax2.spines['right'].set_position(('outward',60))

# ax3 = ax1.twinx()
# ax3.plot(time, Energy, 'blue', aa=True)
#
# ax1.set_ylabel("Thermal Energy Input (Joules)")
# # ax2.set_ylabel("tmperature Injection (Kelvin)")
# ax3.set_ylabel("Heater Energy (Joules)")
# ax1.set_xlabel("Year")
# plt.title("Temperature Constant Injection", fontname = "times new roman")
# plt.grid(True,which='both')
# plt.show()


step1= []
time1 = []
Power1=[]
Energy1=[]
temp1 = []
with open('productionenergy_Kerogen_Temp_Cta.py_700.0.txt', 'r') as file:
    for line in file:
        read = line.split()
        step1.append(float(read[0]))
        time1.append(float(read[1]))
        Power1.append(float(read[2]))
        Energy1.append(float(read[3]))
        temp1.append(float(read[4]))

energy_1 = [] #energy heater 1
energyoini = Energy1[0]
for i in Energy1:
    energy_1.append(energyoini - i)

step2= []
time2 = []
Power2=[]
Energy2=[]
temp2 = []
with open('productionenergy_Kerogen_Temp_Cta.py_800.0.txt', 'r') as file:
    for line in file:
        read = line.split()
        step2.append(float(read[0]))
        time2.append(float(read[1]))
        Power2.append(float(read[2]))
        Energy2.append(float(read[3]))
        temp2.append(float(read[4]))

energy_2 = [] #energy heater 1
energyoini2 = Energy2[0]
for i in Energy2:
    energy_2.append(energyoini2 - i)

step3= []
time3 = []
Power3=[]
Energy3=[]
temp3 = []
with open('productionenergy_Kerogen_Temp_Cta.py_900.0.txt', 'r') as file:
    for line in file:
        read = line.split()
        step3.append(float(read[0]))
        time3.append(float(read[1]))
        Power3.append(float(read[2]))
        Energy3.append(float(read[3]))
        temp3.append(float(read[4]))

energy_3 = [] #energy heater 1
energyoini3 = Energy3[0]
for i in Energy3:
    energy_3.append(energyoini3 - i)


step4= []
time4 = []
Power4=[]
Energy4=[]
temp4 = []
with open('productionenergy_Kerogen_Temp_Cta.py_1000.0.txt', 'r') as file:
    for line in file:
        read = line.split()
        step4.append(float(read[0]))
        time4.append(float(read[1]))
        Power4.append(float(read[2]))
        Energy4.append(float(read[3]))
        temp4.append(float(read[4]))

energy_4 = [] #energy heater 1
energyoini4 = Energy4[0]
for i in Energy4:
    energy_4.append(energyoini4 - i)

#Accumulation Energy

ax = plt.subplot()
ax = plt.plot(time1, temp1, 'blue', linestyle = '-', label= '700 K')
ax1 = plt.plot(time2, temp2, 'green', linestyle = '-', label= '800 K')
ax2 = plt.plot(time3, temp3, 'red', linestyle = '-', label= '900 K')
ax3 = plt.plot(time4, temp4, 'orange', linestyle = '-', label= '1000 K')
plt.grid(True, linewidth=1, which="both")
# plt.legend(loc='upper left')
plt.legend(loc='upper right')
plt.xlabel('Time (Years)')
# plt.ylabel('Energy Input (10^9 Joules)')
plt.ylabel('Injection Power (Joules/sec)')
plt.show()

