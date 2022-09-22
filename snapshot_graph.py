import pygadgetreader as pygr
import matplotlib.pyplot as plt
import numpy as np

p = pygr.readsnap("snap_029", "pos", "dm")
x = np.array(p[:,0])
y = np.array(p[:,1])
print(max(x)-min(x))

p_ic = pygr.readsnap("ic_agora_m12q_ref12_rad6-chull.ics", "pos", "dm")
x_ic = np.array(p_ic[:,0])
y_ic = np.array(p_ic[:,1])
print(max(x)-min(x))

plt.scatter(x, y)
plt.savefig('snap_29_plot.png')
plt.close()

plt.scatter(x_ic, y_ic)
plt.savefig('ic_plot.png')
plt.close()
