import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

f = open("in.fo", "r")

dt = float(f.readline().split(" ")[0])

#crear la figura amb els limits
fig, ax = plt.subplots()
lx = float(f.readline().split(" ")[0])
ly = float(f.readline().split(" ")[0])

x= []
y= []



#llegir el numero de particules i crear contenidors
ntypes = int(f.readline().split(" ")[0])
kE = np.zeros(ntypes)
n = np.zeros(ntypes, int)
for k in range(ntypes):
    n[k] = int(f.readline().split(" ")[0])
    




#settings visuals
col = ['ro', 'bo', 'ko']

while(1):
    ax.clear()
    
    
    #llegir iteracio i parar si hem arribat al final
    it = f.readline()
    if (not it):
        quit()
    
    #llegir energia potencial
    pE = float(f.readline().split(" ")[0])



    #per cada tipus de partícules llegir energia cinetica i posicions
    for k in range(ntypes):
        kE[k] = float(f.readline().split(" ")[0])

        x =[float(i) for i in f.readline().split(" ")]
        y =[float(i) for i in f.readline().split(" ")]
        ax.text(k*130e-10*lx/300e-10,340e-10*ly/300e-10, "KE-"+str(k)+"= "+"{:e}".format(kE[k]) )
        ax.plot(x,y,col[k], markersize=4,)
    
    ax.text(ntypes*130e-10*lx/300e-10,340e-10*ly/300e-10, "PE= "+"{:e}".format(pE) )

    ax.text(20e-10*lx/300e-10,320e-10*ly/300e-10, "temps= "+"{:e}".format(int(it)*dt)+" s   ("+str(int(it))+")" )
    ax.set_xlim([-0.05*lx,1.05*lx])
    ax.set_ylim([-0.05*ly,1.05*ly])
    
    display(fig)
    plt.pause(0.1)




f.close()
