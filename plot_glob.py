import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

plt.close()
data = np.loadtxt('posicions.txt')

Niter = 300;
Npart = 2;
Ntip = [200, 100] 

fig = plt.figure()
a = fig.add_subplot(1,1,1)
lx = 3*(10**-8)
ly = 3*(10**-8)

col = ['ro', 'bo']

for i in range(Niter):
  plt.cla()
  plt.xlim(-0.05*lx, 1.05*lx)
  plt.ylim(-0.05*ly, 1.05*ly)
  counter = 0;
  for j in range(Npart):
      x = data[2*i , counter:counter+Ntip[j]]
      y = data[2*i+1 , counter:counter+Ntip[j]]
      a.plot(x,y,col[j], markersize=4,)
      counter = counter + Ntip[j]
  display(fig)
  plt.pause(0.1)
