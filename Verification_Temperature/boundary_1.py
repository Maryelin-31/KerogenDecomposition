# -*- coding: utf-8 -*-

from zml import *
from zml_hyd import *
from matplotlib import pyplot as plt
import os
import numpy as np

# Create seepage mesh.

mesh = SeepageMesh.create_cube(np.linspace(0, 19.5e-2, 101), np.linspace(0, 13e-2, 26), (-3e-2, 3e-2))

x1, x2 = mesh.get_pos_range(1)
cells_x1 = mesh.get_cells_in_range(xr=(x1 - 0.26e-2, x1 + 0.74e-2))
cells_x2 = mesh.get_cells_in_range(xr=(x2 + 5.76e-2, x2 + 7.76e-2))

y1, y2 = mesh.get_pos_range(1)
cells_y1 = mesh.get_cells_in_range(yr=(y1 - 0.26e-2, y1 + 1.74e-2))
cells_y2 = mesh.get_cells_in_range(yr=(y2 - 1.74e-2, y2 + 0.26e-2))

gas_den = create_ch4_density()
gas_vis = create_ch4_viscosity()

def get_initial_t(x, y, z):
    """
    the initial temperature
    """
    return 274.15 + 24.7

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
    return 0.0, 1.0, 0.0, 0.0

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
model.set(gravity=(0, 0, 0))

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
    cell.set_attr(ca_mc, (vol*2600*1000) - 1.3).set_attr(ca_t, temp).set_attr(ca_g, vol * 1.0)
    cell.fluid_number = 4
    cell.get_fluid(fid_g).set(den=gas_den(pre, temp), vis=gas_vis(pre, temp)).set_attr(fa_t, temp).set_attr(fa_c, 1000)
    cell.get_fluid(fid_w).set(den=1000, vis=1.0e-3).set_attr(fa_t, temp).set_attr(fa_c, 4000)
    cell.get_fluid(fid_o).set(den=oil_den_interp(pre, temp),  vis=oil_vis_interp(pre, temp)).set_attr(fa_t, temp).set_attr(fa_c, 1800)
    cell.get_fluid(fid_k).set(den=2000,  vis=1.0e30).set_attr(fa_t, temp).set_attr(fa_c, 2000).set_attr(ka_dE, 161600.0).set_attr(ka_teq, 565)
    cell.set_attr(ca_k2o, 3600*24*365).set_attr(ca_o2k, 1.0e20)
    cell.fill(pre, sat)
    return cell


for c in mesh.cells:
    cell = add_cell(pos=c.pos, vol=c.vol, pre=get_initial_p(*c.pos), temp=get_initial_t(*c.pos),
                    sat=get_initial_s(*c.pos))
    print(cell.get_attr(ca_g))


for f in mesh.faces:
    face = model.add_face(model.get_cell(f.link[0]), model.get_cell(f.link[1]))
    face.cond = f.area * get_perm(*face.pos) / f.length
    face.set_attr(fa_heatg, f.area * 2.0 / f.length)

tem_env = 24.7 + 274.15 #0 grades Centigrades
# boundary = cells_x1 + cells_x2 + cells_y1 + cells_y2
boundary = cells_y1 + cells_y2
for c in boundary:
    cell = add_cell(pos=c.pos, vol=c.vol, pre=get_initial_p(*c.pos), temp=get_initial_t(*c.pos), sat=get_initial_s(*c.pos))
    cell.set_attr(ca_mc, 1.0).set_attr(ca_t, tem_env)
    face = model.add_face(cell, model.get_cell(c.index))
    face.cond = 0
    face.set_attr(fa_heatg, 3.0)

tem_env = 24.0 + 274.15 #0 grades Centigrades
# boundary = cells_x1 + cells_x2 + cells_y1 + cells_y2
boundary1 = cells_x1 + cells_x2
for c in boundary1:
    cell = add_cell(pos=c.pos, vol=c.vol, pre=get_initial_p(*c.pos), temp=get_initial_t(*c.pos), sat=get_initial_s(*c.pos))
    cell.set_attr(ca_mc, 3.0).set_attr(ca_t, tem_env)
    face = model.add_face(cell, model.get_cell(c.index))
    face.cond = 0
    face.set_attr(fa_heatg, 3.0)


#Injector 1 heater 1
cell2 = model.get_nearest_cell(pos=(9.75e-2, 6.5e-2, 0))
cell2.set_attr(ca_mc, 300)
x, y, z = cell2.pos
rw = 2.5e-2
heat_rate = 10e2 #watts
face2 = model.get_face(cell2.index)
face2.cond = 0
face2.set_attr(fa_heatg, (2 * np.pi * z *2.0) / (np.log(f.length**2) - np.log(rw)))

cell3 = model.get_nearest_cell(pos=(13.85e-2,6.5e-2 , 0))
cell4 = model.get_nearest_cell(pos=(15.45e-2,6.5e-2 , 0))
cell5 = model.get_nearest_cell(pos=(17.05e-2,6.5e-2 , 0))
cell6 = model.get_nearest_cell(pos=(18.65e-2,6.5e-2 , 0))



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

def run_2():
    name = os.path.basename(__file__)
    if os.path.exists(f'production_{name}_{heat_rate}'):
        import shutil
        shutil.rmtree(f'production_{name}_{heat_rate}')
    time = 0
    dt = 1.0
    with open(f'production_{name}_{heat_rate}.txt', 'w') as f:
        for step in range(20000):
            #Temperature of injection
            time += dt
            if time < 2358:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) + (0 * dt) / cell2.get_attr(ca_mc))
            elif time > 2358 and time < 2400:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) + (heat_rate * dt) / cell2.get_attr(ca_mc))
            elif time == 2358 and time == 2400:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) + (heat_rate * dt) / cell2.get_attr(ca_mc))
            elif time > 2400 and time < 3000:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) - (heat_rate * 0.055 * dt) / cell2.get_attr(ca_mc))
            elif time == 2400 or time ==3000:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) - (heat_rate * 0.055  * dt) / cell2.get_attr(ca_mc))
            elif time > 3000 and time < 3480:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) - (heat_rate * 0.009 * dt) / cell2.get_attr(ca_mc))
            elif time == 3000 or time == 3480:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) - (heat_rate * 0.009 * dt) / cell2.get_attr(ca_mc))
            else:
                cell2.set_attr(ca_t, cell2.get_attr(ca_t) - (heat_rate * 0.001 * dt) / cell2.get_attr(ca_mc))

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

            if time > 142*60:
                print(f'stepfinish = {step}')
                break
            if step % 1 == 0:
                print(f'{time/60} {cell3 .get_attr(ca_t)  - 274.15}')
                f.write(f'{time/60} {cell2.get_attr(ca_t) - 274.15} {cell3.get_attr(ca_t) - 274.15} {cell4.get_attr(ca_t) - 274.15} {cell5.get_attr(ca_t) - 274.15} {cell6.get_attr(ca_t) - 274.15}\n')


# run_2()