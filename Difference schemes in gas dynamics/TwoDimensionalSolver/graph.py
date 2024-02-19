import pandas as pd
import os
import sys
import math as m
import matplotlib.pyplot as plt 
from PIL import Image
import numpy as np
from scipy.interpolate import griddata
import time as timeee

def Input(path, **kwarg):
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            for p in list(kwarg.keys()):
                if p == line.split('=')[0].strip():
                    kwarg[p] = float(line.split('=')[1].strip())
    return kwarg.values()


lx, ly, Nx, Ny, variant, cu, tmax = Input('Input.ini', lx=0, ly=0, Nx=0, Ny=0, variant=0, cu=0, tmax=0)

kx = lx/10
ky = ly/1
maxP = 0

print("sys.argv: ", sys.argv)

for i in sys.argv:
    print(i, end="\n")
rank = int(sys.argv[1])
size = int(sys.argv[2])
var = int(variant)
data_path = f'./out/'
time_mass = {}
for adres in os.listdir(data_path):
    if f'Out' not in adres:
        continue
    time_mass[float(adres.split('_')[2])] = adres
    N = int(adres.split('_')[1])
sorted_tuple = sorted(time_mass.items(), key=lambda x: x[0])
time_mass = dict(sorted_tuple)
delta = len(time_mass.keys())%size
step = len(time_mass.keys())//size
if (rank < delta):
    a = rank * step + rank
    b = a + step
else:
    a = rank * step + delta
    b = a + step - 1
if (len(time_mass.keys()) < size):
    if (rank < len(time_mass.keys())):
        a = rank
        b = rank
    else:
        a = 0
        b = -1
p = a
for time in list(time_mass.keys())[a:b + 1]:
    t1 = timeee.time()
    adres = time_mass[time]
    data = pd.read_csv(data_path + adres, sep = ';')
    print(timeee.time() - t1)
    p = data['p']
    rho = data['rho']
    vx = data['vx']
    vy = data['vy']
    x = data['x']
    y = data['y']
    t1 = timeee.time()
    x2 = np.linspace(x.min(),x.max(), m.floor(Nx))
    y2 = np.linspace(y.min(),y.max(), m.floor(Ny))
    xg, yg = np.meshgrid(x2, y2)
    p_z = griddata((x, y), p, (xg, yg), method='linear')
    rho_z = griddata((x, y), rho, (xg, yg), method='linear')
    vx_z = griddata((x, y), vx, (xg, yg), method='linear')
    vy_z = griddata((x, y), vy, (xg, yg), method='linear')
    fig, (axs1, axs2, axs3) = plt.subplots(3, figsize=(10*kx, 4*ky), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    axs1.set_aspect(1)
    axs2.set_aspect(1)
    axs3.set_aspect(1)
    im1 = axs1.pcolormesh(x2, y2, p_z, shading='gouraud', cmap='jet')
    im2 = axs2.pcolormesh(x2, y2, rho_z, shading='gouraud', cmap='jet')
    axs3.streamplot(x2, y2, vx_z, vy_z, density=2.5, linewidth=0.2, arrowstyle='->', arrowsize=0.1, color='r')
    #fig.colorbar(im1, ax=axs1)
    #fig.colorbar(im2, ax=axs2)
    axs1.set_xlim([0, 10])
    axs2.set_xlim([0, 10])
    axs3.set_xlim([0, 10])
    axs1.set_ylim([0, 1])
    axs2.set_ylim([0, 1])
    axs3.set_ylim([0, 1])
    axs1.set_title(f"P(x, y)  t = {time}")
    axs2.set_title(f"Rho(x, y)  t = {time}")
    axs3.set_title(f"V  t = {time}")
    fig.savefig(f"./var{var}/pic/pic_{time:.6f}_.png", dpi = 300)
    p += 1
    fig.clear(False)
    print(timeee.time() - t1)