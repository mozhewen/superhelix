import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pList = np.loadtxt("out.txt")

fig = plt.figure()
ax = Axes3D(fig)
ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([1,1,2,1]))
ax.plot(pList[:,0],pList[:,1],pList[:,2])

plt.show()
