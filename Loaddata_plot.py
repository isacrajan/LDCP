import numpy as np
import matplotlib.pyplot as plt 
import matplotlib 
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

Yc, Uc = np.loadtxt('YcUc.dat', unpack=True)

Xc, Vc = np.loadtxt('XcVc.dat', unpack=True)

plt.figure(1)
plt.plot(Uc,Yc, 'r-', linewidth=2.0)
plt.xlabel(r'\textbf{Centerline U Velocity} (m/s)',fontsize=16)
plt.ylabel(r'\textbf{Y} (m)',fontsize=16)

plt.figure(2)
plt.plot(Xc,Vc, 'b-', linewidth=2.0)
plt.xlabel(r'\textbf{X} (m)',fontsize=16)
plt.ylabel(r'\textbf{Centerline V Velocity} (m/s)',fontsize=16)



plt.show()

