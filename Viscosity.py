
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



#Nuurozieh, 2015 Test
# tempN = [52.7, 62.8, 72.4, 83.0, 92.7, 103.4, 102.9, 127.1, 151.6, 175.5, 190.4] #Nourozieh, 2015 Data Bitumen A
# tempN = [ i + 273.15 for i in tempN]
# tempN = [463.55, 448.65, 424.75, 400.25, 376.05, 376.55, 365.85, 356.15, 345.55, 335.95, 325.85]
# pressure = [2.01, 6.0, 10.0]
# viscosityN = [11.8, 16.2, 32.2, 72.8, 222, 223, 372, 730, 2440, 5560, 14640]
# viscosityN6 = [12.6, 17.5, 34.9, 80.3, 246, 251, 421, 839, 2760, 6700, 17100]
# viscosityN10 =[13.6, 18.8, 37.7, 88.8, 274, 280, 477, 949, 3220, 7710, 19700]
#
# plt.figure()
# plt.plot(tempN,[i * 0.001 for i in viscosityN],'r', label='2.01 MPa')
# plt.plot(tempN,[i * 0.001 for i in viscosityN6],'g', label='6.0 MPa')
# plt.plot(tempN,[i * 0.001 for i in viscosityN10],'b', label='10.0 MPa')
# plt.yscale("log")
# plt.ylabel('viscosity (mPa*sec)')
# plt.xlabel('Temperature (Kelvin)')
# plt.grid(True, linewidth=1, which="both")
# plt.legend()
# plt.show()
# print(tempN)

# Miadonye
# temp = [23.9, 29.9, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 121.0] #Miadonye 1994, Cold Lake bitumen from Mehrotra, 1985b
# temp2 = [i +273.15 for i in temp]
# viscosityE = [54.50, 24.10, 7.64, 2.97, 1.34, 0.66, 0.35, 0.21, 0.13, 0.09, 0.06]
# viscosityM = [52.83, 24.99, 8.17, 3.14, 1.37, 0.66, 0.35, 0.20, 0.12, 0.08, 0.05]
# plt.plot(temp2,viscosityE, 'r', label='Experimental Viscosity')
# plt.plot(temp2,viscosityM, 'g', label='Predicted Viscosity')
# plt.yscale("log")
# plt.ylabel('viscosity (Pa*sec)')
# plt.xlabel('Temperature (Kelvin)')
# plt.xlim(280,400)
# plt.grid(True, linewidth=1, which="both")
# plt.legend()
# plt.show()


# temp = np.linspace(300, 600, 100)
# # pressure = np.linspace(1.0, 40.0, 100)
# pressure = 0.1 #MPa

# def oil_vis(temp): #Miadonye, 1994
#     c = -3.0020
#     b = np.log(29.9) - c #viscosity at 30C == 29.9 Pa*sec at 1 atm
#     s = 0.0066940 * b + 3.5364
#     Vis = []
#     for i in temp:
#         vis = np.exp((b /((1 + ((i - 273.15) - 30)/303.15)**s)) + c)
#         Vis.append(vis)
#     return Vis
#
# def oil_vis(temp): #Miadonye, 1994
#     c = -3.0020
#     b = np.log(24) - c #viscosity at 30C == 29.9 Pa*sec at 1 atm
#     s = 0.0066940 * b + 3.5364
#     vis = np.exp((b /((1 + ((temp- 273.15) - 30)/303.15)**s)) + c)
#     return vis

#
# # def oilVis(pressure, temp): #Mehrotra and Svrcek, 1986
# #     b1=24.84525
# #     b2=-3.90450
# #     b3=0.004723
# #     for i in pressure:
# #         A = (b1 + (b2 * np.log(temp))) + (b3 * i)
# #         vis_oil = 0.001 * (np.exp(np.exp(A))) #Pa*sec
# #     return vis_oil
#
# def oil_Vis(pressure, temp): #Mehrotra and Svrcek, 1986
#     b1=22.8515
#     b2=-3.5784
#     b3=0.00511938
#     A = (b1 + (b2 * np.log(temp))) + (b3 *(pressure))
#     vis_oil = 0.001 * (np.exp(np.exp(A))) #Pa*sec
#     return vis_oil
#
# print(oil_vis(356))
# print(oil_Vis(3.5,356))

# plt.figure()
# plt.plot(temp,oil_vis(temp), 'g', label='Miadonye Model')
# plt.plot(temp,oil_Vis(pressure, temp),'r', label='Mehrotra Model')
# plt.yscale("log")
# plt.ylabel('viscosity (Pa*sec)')
# plt.xlabel('Temperature')
# plt.grid(True, linewidth=1, which="both")
# plt.legend()
# plt.show()


# DENSITY


#
# def oil_den(pressure, temp):
#     a1 = 1021.62 #kg*m^-3
#     a2 = - 0.58976 #kg*m^-3 * C^-1
#     a4 = 0.382 #1/MPa
#     a5 = 0.00283 #C^-1
#     alpha = a4 + a5 * (temp - 273.5) #temp = Kelvin to Celsius
#     den_o = a1 + a2 * (temp - 273.15)
#     den = den_o + alpha * (pressure * 0.000001)#Pressure Mpa
#     return den

# print(oil_den(pressure, temp))







