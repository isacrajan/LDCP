import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

plt.style.use('seaborn-deep')
#List all available styles: print(plt.style.available)
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 

Kn = np.array([0.01,0.05,0.1,0.15])
#sig = np.array([0.2, 0.5, 0.7, 1])
mark = ['-', ':', '--', '-.']

fig = plt.figure(1)
for i in range(len(Kn)):
    sig = 0.7
    a = mark[i]
    name = "UcYc_Kn" + repr(round(Kn[i],2)) + "sig" + repr(round(sig,1)) + ".plt"
    y, u = np.loadtxt(name, unpack=True)
    plt.tick_params(labelsize=13, direction='in')
    plt.plot(y,u, mark[i],linewidth=1.5, label='Kn ' + repr(round(Kn[i],3)) + ' TMAC ' + repr(round(sig,1)) )
    
plt.xlabel('Centerline u/U',fontsize=15)
plt.ylabel('Y (m)',fontsize=15)
plt.axis([-0.2, 1.1, -0.1, 2.1])
#plt.axis('tight')
plt.grid()
#plt.title("Shallow Cavity, K=2")
plt.legend(loc='best', prop={'size':11}, frameon=True, shadow=True, framealpha=0.5) #bbox_to_anchor=(1.1, 1.1)
plt.show()
fig.savefig('Y vs U sig ' + repr(round(sig,1)) + '.pdf')    

fig = plt.figure(2)
for i in range(len(Kn)):
    sig = 0.7
    a = mark[i]
    name = "XcVc_Kn" + repr(round(Kn[i],2)) + "sig" + repr(round(sig,1)) + ".plt"
    x, v = np.loadtxt(name, unpack=True)
    plt.tick_params(labelsize=13, direction='in')
    plt.plot(x,v, mark[i],linewidth=1.5, label='Kn ' + repr(round(Kn[i],3)) + ' TMAC ' + repr(round(sig,1)) ) 
plt.xlabel('X (m)',fontsize=15)
plt.ylabel('Centerline v/U',fontsize=15)
plt.legend(loc='best',prop={'size':11}, frameon=True, shadow=True, framealpha=0.5)
plt.axis([-0.1, 1.1, -0.5, 0.5])
#plt.axis('tight')
plt.grid()
#plt.title("Shallow Cavity, K=2")
plt.show()
fig.savefig('V vs X sig ' + repr(round(sig,1)) + '.pdf')

#plt.close("all")