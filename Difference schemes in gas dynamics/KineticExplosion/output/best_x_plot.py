import matplotlib.pyplot as plt
import pandas as pd
import os

# Plot of le and best_x depends
# read data from result.csv
df = pd.read_csv('result.csv', sep = ';')

le = df['le']
best_x = df['best_x']

# plot
plt.scatter(le, best_x, color='blue', marker='o', s=5)
plt.title('Plot of le and zone, where pressure is bigger than 1.5atm')
plt.xlabel('le')
plt.ylabel('best_x')
plt.savefig('best_x_plot.png')
plt.show()