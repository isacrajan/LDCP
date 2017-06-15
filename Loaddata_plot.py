import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-deep')
#List all available styles: print(plt.style.available)
plt.rc('font', family='serif') 
plt.rc('font', serif='Palatino Linotype') 

#Kn 0.01 and sigma 0.7
Yc1, Uc1 = np.loadtxt('YcUc_Kn0.01sig0.7.dat', unpack=True)
Xc1, Vc1 = np.loadtxt('XcVc_Kn0.01sig0.7.dat', unpack=True)
#Kn 0.05 and sigma 0.7
Yc3, Uc3 = np.loadtxt('YcUc_Kn0.05sig0.7.dat', unpack=True)
Xc3, Vc3 = np.loadtxt('XcVc_Kn0.05sig0.7.dat', unpack=True)
#Kn 0.1 and sigma 0.7
Yc5, Uc5 = np.loadtxt('YcUc_Kn0.1sig0.7.dat', unpack=True)
Xc5, Vc5 = np.loadtxt('XcVc_Kn0.1sig0.7.dat', unpack=True)

#U_ben, Y_ben = np.loadtxt('Kn0.01sig0.7 YvsU.dat', unpack=True)
#X_ben, V_ben = np.loadtxt('Kn0.01sig0.7 VvxX.dat', unpack=True)

fig = plt.figure(1)
plt.tick_params(labelsize=12)
plt.plot(Uc1,Yc1, '-', linewidth=2.0, label='Kn 0.01 TMAC 0.7')
plt.plot(Uc3,Yc3, '--', linewidth=2.0, label='Kn 0.05 TMAC 0.7')
plt.plot(Uc5,Yc5, ':', linewidth=2.0, label='Kn 0.1 TMAC 0.7')
#plt.scatter(U_ben, Y_ben, marker='o', s=35, alpha=0.5)
plt.xlabel('Centerline u/U',fontsize=15)
plt.ylabel('Y (m)',fontsize=15)
plt.legend(loc='best',prop={'size':11}, frameon=True)
plt.grid()
plt.axis([-0.2, 1, -0.1, 1.1])
fig.savefig('Y vs U.pdf')

fig = plt.figure(2)
plt.tick_params(labelsize=12)
plt.plot(Xc1,Vc1, '-', linewidth=2.0, label='Kn 0.01 TMAC 0.7')
plt.plot(Xc3,Vc3, '--', linewidth=2.0, label='Kn 0.05 TMAC 0.7')
plt.plot(Xc5,Vc5, ':', linewidth=2.0, label='Kn 0.1 TMAC 0.7')
#plt.scatter(X_ben, V_ben, marker='o', s=35, alpha=0.5)
plt.xlabel('X (m)',fontsize=15)
plt.ylabel('Centerline v/U',fontsize=15)
plt.legend(loc='best',prop={'size':11}, frameon=True)
plt.grid()
plt.axis([-0.1, 1.1, -0.15, 0.15])
fig.savefig('U vs X.pdf')

plt.show()

#plt.close("all")