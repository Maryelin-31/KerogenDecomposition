# -*- coding: utf-8 -*-

from zml import *
from zml_hyd import *
from matplotlib import pyplot as plt
from matplotlib import ticker, cm
from matplotlib.ticker import LinearLocator
from matplotlib.cm import ScalarMappable
import os
import numpy as np




# Create seepage mesh.

mesh = SeepageMesh.create_cube(np.linspace(0, 50, 101), np.linspace(0, 50, 26), (-1.5, 1.5))


gas_den = create_ch4_density()
gas_vis = create_ch4_viscosity()

def get_initial_t(x, y, z):
    """
    the initial temperature
    """
    return 338.0 + 22.15 - 0.0443 * y

def get_initial_p(x, y, z):
    """
    the initial pressure
    """
    return 15.0e6 + 5e6 - 1e4 * y

def get_perm(x, y, z):
    """
    the initial permeability
    """
    return 1.0e-15

def get_initial_s(x, y, z):
    """
    the initial saturation (gas, water, oil, kerogen)
    """
    return 0.1, 0.0, 0.3, 0.6

# The ID of gas, water, oil and kerogen
fid_g = 0
fid_w = 1
fid_o = 2
fid_k = 3

# The ID of temperature and specific heat
fa_t = 0
fa_c = 1

# kerogen attrbution:
#  ka_dE: the energy needed for convert 1kg kerogen to oil
#  ka_teq: the temperature above which kerogen can be converted to oil
ka_dE = 2
ka_teq = 3

# Cell attributions
ca_vol = 0   # Cell volume
ca_mc = 1    # rock mass multiply heat capacity in a cell
ca_t = 2     # rock temperature
ca_fp = 3    # cell fluid pressure
ca_g = 4     # the conductivity for exchange heat between rock and fluids
ca_k2o = 5   # the half time for convert kerogen to oil
ca_o2k = 6   # the half time for convert oil to kerogen (set it to infinite value to disable kerogen to oil)

# Face Attribution
fa_heatg = 0 # the heat 0000000 for face.

model = Seepage()
model.set(gravity=(0, -10, 0))


def oil_vis(pressure, temp): #Mehrotra and Svrcek, 1986
    b1 = 22.8515
    b2 = -3.5784
    b3 = int(0.00511938)
    A = (b1 + (b2 * np.log(temp))) + (b3 * (pressure * 0.000001))
    vis_oil = 0.001 * (np.exp(np.exp(A)))
    return vis_oil

def create_oil_viscosity():
    data = Interp2()
    data.create(1.0e6, 0.1e6, 40e6, 300, 1, 1200, oil_vis)
    return data

def oil_den(pressure, temp): 
    a1 = 1021.62 #kg*m^-3
    a2 = - 0.58976 #kg*m^-3 * C^-1
    a4 = 0.382 #1/MPa
    a5 = 0.00283 #C^-1
    alpha = a4 + a5 * (temp - 273.5) #temp = Kelvin to Celsius
    den_o = a1 + a2 * (temp - 273.15)
    den = den_o + alpha * (pressure * 0.000001)#Pressure Mpa
    return den

def create_oil_density():
    data = Interp2()
    data.create(1.0e6, 0.1e6, 40e6, 300, 1, 1200, oil_den)
    return data

oil_den_interp = create_oil_density()
oil_vis_interp = create_oil_viscosity()

def add_cell(pos, vol, pre, temp, sat):
    cell = model.add_cell().set(pos=pos).set_attr(ca_vol, vol).set_pore(10e6, vol * 0.43, 1000.0e6, vol)
    cell.set_attr(ca_mc, vol*2600*1000).set_attr(ca_t, temp).set_attr(ca_g, 100) #heat transfer coefficient
    cell.fluid_number = 4
    cell.get_fluid(fid_g).set(den=gas_den(pre, temp), vis=gas_vis(pre, temp)).set_attr(fa_t, temp).set_attr(fa_c, 2000)
    cell.get_fluid(fid_w).set(den=1000, vis=1.0e-3).set_attr(fa_t, temp).set_attr(fa_c, 4200)
    cell.get_fluid(fid_o).set(den=oil_den_interp(pre, temp),  vis=oil_vis_interp(pre, temp)).set_attr(fa_t, temp).set_attr(fa_c, 1800)
    cell.get_fluid(fid_k).set(den=1500,  vis=1.0e30).set_attr(fa_t, temp).set_attr(fa_c, 2000).set_attr(ka_dE, 161600.0).set_attr(ka_teq, 565)
    cell.set_attr(ca_k2o, 3600*24*(365*1)).set_attr(ca_o2k, 1.0e20)
    cell.fill(pre, sat)
    return cell


for c in mesh.cells:
    cell = add_cell(pos=c.pos, vol=c.vol, pre=get_initial_p(*c.pos), temp=get_initial_t(*c.pos),
                    sat=get_initial_s(*c.pos))


for f in mesh.faces:
    face = model.add_face(model.get_cell(f.link[0]), model.get_cell(f.link[1]))
    face.cond = f.area * get_perm(*face.pos) / f.length
    face.set_attr(fa_heatg, f.area * 1.0 / f.length)
    

#Producer
cell_1 = model.get_nearest_cell(pos=(25, 25, 0))
vol = 1.0e6 #increase the cell volume
rw = 3.0 #radious o well
cell1 = add_cell(pos=(25,25,1), vol=vol, pre=get_initial_p(25,25,1), temp=get_initial_t(25,25,1), sat=get_initial_s(25,25,1))
x, y, z = cell1.pos #declare the width z
cell1.set_pore(10.0e6, vol * 0.1, 1000.0e6, vol * 0.1)
cell1.fill(3.5e6, get_initial_s(*cell1.pos))
cell1.set_attr(ca_g, 1.0e-20)
face = model.add_face(cell1, cell_1)
face.cond = (2 * np.pi * z *  get_perm(*face.pos)) / (np.log(f.length**2) - np.log(rw)) #radial flow
face.set_attr(fa_heatg, 0)

#Monitor Production Well
def monitor(fid):
    mass = cell1.get_fluid(fid).mass      
    volume = cell1.get_fluid(fid).vol
    saturation = cell1.get_fluid(fid).vol_fraction
    return mass, volume, saturation
#oil
mass_ini, volume_ini, saturation_ini = monitor(fid_o)

#Injector 1 heater 1
cell2 = model.get_nearest_cell(pos=(20, 25, 0))
heat_rate = 2000 #watts 
face2 = model.get_face(cell2.index)
face2.cond = 0
face2.set_attr(fa_heatg, (2 * np.pi * z * 1.0) / (np.log(f.length**2) - np.log(rw)))

#Injector 2 heater 2
cell3 = model.get_nearest_cell(pos=(30, 25, 0))
face3 = model.get_face(cell3.index)
face3.cond = 0
face3.set_attr(fa_heatg, (2 * np.pi * z * 1.0) / (np.log(f.length**2) - np.log(rw)))

def total_mass(fid):
    mass = []
    for cell in model.cells:
        mass.append(cell.get_fluid(fid).mass)
    mass = mass[:len(mass)-1]
    total_mass = sum(mass)
    return total_mass

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
    return sw, krw, sg, krg, so, kro
sw, krw, sg, krg, so, kro = stone_model_I(swir=0.1, sorg=0.1, sorw=0.1, sgc=0.1, krwro=0.9, kroiw=1, krgro=0.9, nw=2, nsorw=2, ng=2, nog=2)

# Set relative permeability.
model.set_kr(fid_g, sg, krg)
model.set_kr(fid_w, sw, krw)
model.set_kr(fid_o, so, kro)

def saturation(fid, time, step):
    name = os.path.basename(__file__)
    X = []
    Y = []
    Saturation = []
    for cell in model.cells:
        x, y, z = cell.pos
        X.append(x)
        Y.append(y)
        Cell_Fluid_Vol = 0
        for fluid in range(cell.fluid_number):
            Cell_Fluid_Vol += cell.get_fluid(fluid).vol
        Cell_Fluid = cell.get_fluid(fid).vol
        saturation = Cell_Fluid / Cell_Fluid_Vol
        Saturation.append(saturation)
    saturation = Saturation[: len(Saturation)-1]
    Sat = np.array(saturation)
    S = np.transpose(Sat.reshape((Sat.size//25,25)))
    [xx, yy] = np.meshgrid(X, Y)
    fig, ax = plt.subplots()
    if fid == fid_g:
        fid = 'gas'
        vmin = get_initial_s(x, y, z)[0]
        vmax= 1.0
    if fid == fid_w:
        fid = 'water'
        sat_ini = get_initial_s(x, y, z)[1]   
    if fid == fid_o:
        fid = 'oil'
        vmin = get_initial_s(x, y, z)[2]
        vmax= 1.0
    if fid == fid_k:
        fid = 'kerogen'
        vmin = 0
        vmax=get_initial_s(x, y, z)[3]
    plot = plt.contourf(S, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap='summer', vmin = vmin, vmax = vmax)
    Time = round(time / (3600*24*365), 1)
    plt.title(f'Saturation of {fid} at {Time} year ', fontsize=10)
    plt.xlabel('x/m')
    plt.ylabel('y/m')
    plt.clim(vmin=vmin, vmax=vmax)
    bar = np.linspace(vmin, vmax, int(vmax * 10))
    fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap), ticks=bar)

    PlotPath = f'Saturation_fluid_{fid}_step_{step}.png'
    folder = os.path.join(os.getcwd(), f'Saturation_fluid_Plot_{name}')
    plotfile = os.path.join(folder, PlotPath)
    plt.savefig(plotfile , format='png')
    # plt.show()
    plt.close()

def Temp_cell(time,step): #Matrix 
    name = os.path.basename(__file__)
    X = []
    Y = []
    Temp = []
    for cell in model.cells:
        x, y, z = cell.pos
        X.append(x)
        Y.append(y)
        Temp.append(cell.get_attr(ca_t))  
    temp = Temp[:len(Temp) - 1]
    [xx, yy] = np.meshgrid(X, Y)
    fig, ax = plt.subplots()
    temperature = np.array(temp)
    Temperature= np.transpose(temperature.reshape((temperature.size // 25, 25)))
    tini = 300
    tmax = np.max(temperature)
    plot = plt.contourf(Temperature, 20, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap='coolwarm', vmin=tini, vmax=tmax)
    Time = round(time / (3600*24*365), 1)
    plt.title(f'Temperature of Rock (K) at {Time} year', fontsize=10)
    plt.xlabel('x/m')
    plt.ylabel('y/m')
    plt.clim(vmin=tini, vmax=tmax)
    fig.colorbar(ScalarMappable(norm=plot.norm, cmap=plot.cmap))
    
    PlotPath = f'Temperature_{step}.png'   
    folder = os.path.join(os.getcwd(),f'Temp_cell_Plot_{name}')
    plotfile = os.path.join(folder, PlotPath)
    plt.savefig(plotfile , format='png')
    # plt.show()
    plt.close()


def Save(path, step, time):
    name = os.path.basename(__file__)
    assert isinstance(path,str)
    Savefolder = os.path.join(os.getcwd(), f'Model_Results_{name}_{heat_rate}')
    SavePath = os.path.join(Savefolder, path)
    with open(SavePath, 'w') as file:
        Time = round(time / (3600*24*365), 2)
        for cell in model.cells:
            x, y, z = cell.pos
            file.write(f'{step} {Time} {x} {cell.get_fluid(fid_k).mass} {cell.get_fluid(fid_o).mass} {cell.get_fluid(fid_k).get_attr(fa_t)} {cell.get_attr(ca_t)} {cell.get_fluid(fid_k).get_attr(ka_teq)}\n')
             
def run_2():    
    name = os.path.basename(__file__)
    if os.path.exists(f'Model_Results_{name}_{heat_rate}'):
        import shutil
        shutil.rmtree(f'Model_Results_{name}_{heat_rate}')
    os.mkdir(f'Model_Results_{name}_{heat_rate}')
    if os.path.exists(f'production_{name}_{heat_rate}'):
        import shutil
        shutil.rmtree(f'production_{name}_{heat_rate}')
    # if os.path.exists(f'Saturation_fluid_Plot_{name}'):
    #     import shutil
    #     shutil.rmtree(f'Saturation_fluid_Plot_{name}')
    # os.mkdir(f'Saturation_fluid_Plot_{name}')
    # if os.path.exists(f'Temp_cell_Plot_{name}'):
    #     import shutil
    #     shutil.rmtree(f'Temp_cell_Plot_{name}')
    # os.mkdir(f'Temp_cell_Plot_{name}')
    time = 0
    dt = 1.0e-8  
    # cell2ini = cell2.get_attr(ca_t)
    # cell3ini = cell3.get_attr(ca_t)
    with open(f'production_{name}_{heat_rate}.txt', 'w') as f:
        for step in range(20000):
            #Temperature of injection
            cell2.set_attr(ca_t, cell2.get_attr(ca_t) + (heat_rate * dt) / cell2.get_attr(ca_mc))    
            cell3.set_attr(ca_t, cell3.get_attr(ca_t) + (heat_rate * dt) / cell3.get_attr(ca_mc))
            model.iterate(dt, ca_p=ca_fp)
            model.iterate_thermal(ca_t=ca_t, ca_mc=ca_mc, fa_g=fa_heatg, dt=dt)
            model.exchange_heat(dt=dt, ca_g=ca_g, ca_t=ca_t, ca_mc=ca_mc, fa_t=fa_t, fa_c=fa_c)
            
            update_ice(model, dt=dt, fid_i=fid_k, fid_w=fid_o, fa_t=fa_t, fa_c=fa_c, ia_dE=ka_dE, ia_teq=ka_teq,
                        ca_i2w=ca_k2o, ca_w2i=ca_o2k)          
            model.update_den(fluid_id=fid_g, kernel=gas_den, relax_factor=0.01,
                              fa_t=fa_t, min=1, max=100)
            model.update_den(fluid_id=fid_o, kernel=oil_den_interp, relax_factor=0.01,
                              fa_t=fa_t, min=800, max=1000)            
            if step % 1 == 0:
                model.update_vis(fluid_id=fid_g, kernel=gas_vis, ca_p=ca_fp, fa_t=fa_t,
                                 relax_factor=0.1, min=1.0e-6, max=1.0e-3)

                model.update_vis(fluid_id=fid_o, kernel=oil_vis_interp, ca_p=ca_fp, fa_t=fa_t,
                                 relax_factor=0.1, min=1.0e-2, max=20)
            time += dt
            if time > 3600*24*365*20:
                print(f'stepfinish = {step}')
                break
            dt = model.get_recommended_dt(dt)
            dt = max(1.0e-6, min(3600*24*5, dt))
            path = f'step{step}.txt' 
            if step % 1 == 0:
                Save(path, step, time)
                # saturation(fid_k, time, step)
                # Temp_cell(time,step)
                print(f'{round(time/(3600*24*365),1)}, {step}')
                f.write(f'{step} {round(time/(3600*24*365),1)} {cell1.get_fluid(fid_o).mass -  mass_ini} {total_mass(fid_k)} {total_mass(fid_o)} {cell2.get_attr(ca_t)} {cell3.get_attr(ca_t)} {cell_1.get_fluid(fid_o).get_attr(fa_t)}\n')
                
                
run_2()

