import matplotlib.pyplot as plt
import pandas as pd
import os
from PIL import Image


def find_max_min(data, s):
    values = []
    for k in data.keys():
        for value in data[k][s]:
            try:
                numeric_value = float(value)
                values.append(numeric_value)
            except ValueError:
                print(f"Ignoring non-numeric value: {value}")

        # l += list(data[k][s])

    return max(values), min(values)

path = os.path.join(os.getcwd(), 'out', 'csv')
data = {}
for file in os.listdir(path):
    t = file.split('-')[1]
    data[float(t)] = pd.read_csv(path + f'/{file}', sep = ';')

p_max, p_min = find_max_min(data, 'p')
v_max, v_min = find_max_min(data, 'v')
r_max, r_min = find_max_min(data, 'r')
e_max, e_min = find_max_min(data, 'e')


for k in data.keys():

    p = data[k]['p']
    v = data[k]['v']
    r = data[k]['r']
    e = data[k]['e']
    x = data[k]['x']
    #define grid of plots
    #add title
    fig, axs = plt.subplots(nrows= 2 , ncols= 2)
    fig.suptitle(f'Time = {k}')
    #add data to plots

    axs[0, 0].plot(x, p, 'ko', markersize=1)
    axs[0, 1].plot(x, v, 'ko', markersize=1)
    axs[1, 0].plot(x, r, 'ko', markersize=1)
    axs[1, 1].plot(x, e, 'ko', markersize=1)

    axs[0, 0].set_title('p')
    axs[0, 0].set_ylim([p_min, p_max])
    axs[0, 1].set_title('v')
    axs[0, 1].set_ylim([v_min, v_max])
    axs[1, 0].set_title('r')
    axs[1, 0].set_ylim([r_min, r_max])
    axs[1, 1].set_title('e')
    axs[1, 1].set_ylim([e_min, e_max])

    fig.tight_layout()
    fig.savefig(f'./out/pic/pic_{k}_.jpg', dpi = 150)
    fig.clf()

sorted_dict = dict(sorted(data.items()))
duration = 15000/len(sorted_dict)
image_path_list = [f'./out/pic/pic_{k}_.jpg' for k in sorted_dict.keys()]
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