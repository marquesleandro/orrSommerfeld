import numpy as np
import matplotlib.pyplot as plt
from math import pi
import pylab as pl

#eq de blausius obtida a partir da espessura da camada limite atraves da solucao do problema de Rayleigh

deta=0.002
T=10
Pr=1        #1      100.    0.1

#deta=0.005
#T=    1          2       3         10
#h[0]=1.015      0.6     0.48      0.46 



eta = np.linspace(0, T, int(T/deta)+1) # An array to store the solution
f = np.zeros(len(eta))
g = np.zeros(len(eta))
h = np.zeros(len(eta))

def df(g):
			return g

def dg(h):
			return h
			
def dh(f,h):
			return -f*h
			
f[0] = 0
g[0] = 0
# GUESS h(1) SUCH THAT g(N)=1, f'(N)=1
#quando mudar o tempo o ideal seria ajustar o h automaticamente...mathematica ok
h[0]=0.46     
			
for i in range(1,len(eta)):
		f[i] = f[i-1] + df(g[i-1])*deta
		g[i] = g[i-1] + dg(h[i-1])*deta
		h[i] = h[i-1] + dh(f[i-1],h[i-1])*deta
		
		#if abs(g[i] - 1)<=10**(-6):
        #	print(f[i])
        #	break
		
#np.savetxt('eta.txt', eta)  #Save the solution
#np.savetxt('f.txt', f)
#np.savetxt('df.txt', g)
#np.savetxt('d2f.txt', h)
plt.figure(figsize=(12,6))

#plt.subplot(211)
plt.subplot(221)
plt.title(r'Blausius equation (boundary layer)', fontsize=15)
plt.plot(f,eta,color='blue')
#plt.plot(eta, h,color='purple')
plt.ylabel('$\eta$ axis')
plt.xlabel('f axis ')
plt.legend(('f', 'g=df'), loc='upper left') #, 'h=d2f'
plt.text(3., 3,'Blausius equation: d3f+f*d2f=0', fontsize=12),
plt.grid(True)

#plt.subplot(212)
plt.subplot(223)
plt.text(2, 0.67, r'objetivo:', fontsize=9)
plt.text(2, 0.6, r'achar f com a c.c df(N)=1. Definindo $\eta$=10, atraves', fontsize=9)
plt.text(2, 0.53, r'do metodo do chute, quando d2f(0)=0.46, df(N)=1.', fontsize=9),
plt.plot(eta, g,color='green')
#plt.plot(eta,f,'o',color='blue')
#plt.plot(eta, h,color='purple')
plt.xlabel('$\eta$, $\eta$:0:10, d$\eta$=0.002')
plt.ylabel('df axis ')
plt.legend(('df','f'), loc='upper left') #, 'h=d2f'
plt.grid(True)

plt.subplot(222)
plt.title(r'Blasius flow profile', fontsize=15)
plt.plot(g,eta,color='green')
plt.plot(f,eta,color='blue')
plt.text(2.5, 9,'df=Vx/U', fontsize=12),
plt.ylabel('$\eta$, $\eta$:0:10, d$\eta$=0.002')
#plt.xlabel('df=Vx/U axis ')
plt.legend(('df','f'), loc='upper left') 
plt.grid(True)
	
plt.savefig('blausius-1.pdf')
plt.show()
