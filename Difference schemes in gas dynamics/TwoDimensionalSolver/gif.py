import pandas as pd
import os
import sys
import math as m
import matplotlib.pyplot as plt 
from PIL import Image
import numpy as np

def Input(path, **kwarg):
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            for p in list(kwarg.keys()):
                if p == line.split('=')[0].strip():
                    kwarg[p] = float(line.split('=')[1].strip())
    return kwarg.values()
lx, ly, Nx, Ny, variant, cu, tmax = Input('Input.ini', lx=0, ly=0, Nx=0, Ny=0, variant=0, cu=0, tmax=0)

var = int(variant)
data_path = f'./var{var}/pic/'
time_mass = {}
for adres in os.listdir(data_path):
    if f'pic' not in adres:
        continue
    time_mass[float(adres.split('_')[1])] = adres
sorted_tuple = sorted(time_mass.items(), key=lambda x: x[0])
time_mass = dict(sorted_tuple)

duration_l = list(time_mass.keys())
frames = []
duration_l = (np.array(duration_l[1:]) - np.array(duration_l[:-1]))
duration_l = np.append(duration_l, duration_l[len(duration_l) - 1])
ctime = 15 
tscale = None
p = 0
for time in list(time_mass.keys()):
    frame = Image.open(data_path + time_mass[time])
    p += 1
    frames.append(frame)
if tscale == None:
    tscale = ctime * 1E3 / duration_l.sum()
    duration_l = duration_l * tscale
duration_list = []
for i in duration_l:
    i = m.floor(i)
    duration_list.append(i)
frames[0].save(
f'./var{var}/anime.gif',
save_all=True,
append_images=frames[1:],
optimize=True,
duration=duration_list,
loop=0
)