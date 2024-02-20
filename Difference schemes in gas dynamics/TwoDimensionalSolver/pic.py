from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from PIL import Image

path = './out/csv'
data = {}
for file in os.listdir(path):
    t = file.split('_')[1]
    data[float(t)] = pd.read_csv(path + f'/{file}', sep = ';')

data = dict(sorted(data.items()))
print(len(data))

for k in data.keys():
    x = data[k]['x']
    y = data[k]['y']
    p = data[k]['p']
    vx = data[k]['vx']
    vy = data[k]['vy']
    r = data[k]['r']
    e = data[k]['e']

    len_x = len(set(x))
    len_y = len(set(y))
    print(len_x, len_y)
    x_mash = np.linspace(x.min(),x.max(), len_x)
    y_mash = np.linspace(y.min(),y.max(), len_y)
    xg, yg = np.meshgrid(x_mash, y_mash)

    p_mesh = griddata((x, y), p, (xg, yg), method='linear')
    vx_mesh = griddata((x, y), vx, (xg, yg), method='linear')
    vy_mesh = griddata((x, y), vy, (xg, yg), method='linear')
    r_mesh = griddata((x, y), r, (xg, yg), method='linear')

    fig, axs = plt.subplots(nrows= 4 , ncols= 1)
    fig.suptitle(f'Time = {k}')

    axs[0].pcolormesh(x_mash, y_mash, p_mesh, shading='gouraud', cmap='jet')
    axs[0].set_aspect('equal', adjustable='box')
    axs[0].set_title('p')

    axs[1].pcolormesh(x_mash, y_mash, vx_mesh, shading='gouraud', cmap='jet')
    axs[1].set_aspect('equal', adjustable='box')
    axs[1].set_title('vx')

    axs[2].pcolormesh(x_mash, y_mash, vy_mesh, shading='gouraud', cmap='jet')
    axs[2].set_aspect('equal', adjustable='box')
    axs[2].set_title('vy')

    axs[3].pcolormesh(x_mash, y_mash, r_mesh, shading='gouraud', cmap='jet')
    axs[3].set_aspect('equal', adjustable='box')
    axs[3].set_title('r')

    fig.tight_layout()
    fig.savefig(f'./out/pic/pic_{k}_.jpg', dpi = 150)
    fig.clf()

duration = 15000/len(data)
image_path_list = [f'./out/pic/pic_{k}_.jpg' for k in data.keys()]
image_list = [Image.open(file) for file in image_path_list]
image_list[0].save(
            r'./out/anime.gif',
            save_all=True,
            append_images=image_list[1:], # append rest of the images
            duration=duration, # in milliseconds
            loop=0)

i = 1
while (True):
    if not os.path.exists(f'./out({i})'):
        os.mkdir(f'./out({i})')
        break
    else:
        i += 1

new_path = f'./out({i})'
os.rename('./out/csv', new_path + '/csv')
os.rename('./out/pic', new_path + '/pic')
os.rename('./out/anime.gif', new_path + '/anime.gif')

os.mkdir('./out/csv')
os.mkdir('./out/pic')