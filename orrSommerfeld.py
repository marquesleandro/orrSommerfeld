import numpy as np
from numpy import linalg
from scipy.linalg import eigvals
import matplotlib.pyplot as plt
from math import pi
import pylab as pl



# -----------------
# BLAUSIUS EQUATION
# -----------------

# --------------------------------------------------------------------------------------------
deta=0.002
T=6

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
			
f[0] = 0.0
g[0] = 0.0
h[0] = 0.47   
			
for i in range(1,len(eta)):
		f[i] = f[i-1] + df(g[i-1])*deta
		g[i] = g[i-1] + dg(h[i-1])*deta
		h[i] = h[i-1] + dh(f[i-1],h[i-1])*deta
		
plt.title(r'Blasius flow profile', fontsize=15)
plt.plot(eta,g,color='black')
plt.xlabel('$\eta$')
plt.ylabel('$df$')
plt.grid(True)
	
#plt.savefig('blausius-1.pdf')
plt.show()
# --------------------------------------------------------------------------------------------







# -----------------------
# ORR-SOMMERFELD EQUATION
# -----------------------


# --------------------------------------------------------------------------------------------
numElements = 50
numNodes = numElements + 1
L = 1.0
y = np.zeros([numNodes,1], dtype=float)
yDelta = L/numElements
for i in range(0,numNodes):
 y[i] = (L/numElements)*i


u = np.zeros([numNodes,1], dtype=float)
u2 = np.zeros([numNodes,1], dtype=float)
uMax = 1.5
for i in range(0,numNodes):
 u[i] = 4.0*uMax*((y[i]/L) - ((y[i]/L)**2))
 u2[i] = -8.0*uMax/(L**2)


Re = 12000
alpha = np.linspace(1.0,2.4,num=10)
c = np.zeros([len(alpha),1],dtype=float)

for k in range(0,len(alpha)):

 # Christopher Simpson: 
 # https://www.cdsimpson.net/2015/04/temporal-and-spatial-stability-analysis.html
 A = np.zeros([numNodes,numNodes], dtype=complex)
 a1 = (-u[1]*alpha[k]**3) - (u2[1]*alpha[k]) + ((1j*alpha[k]**4)/Re)
 a2 = ( u[1]*alpha[k]) - ((2.0*1j*alpha[k]**2)/Re)
 a3 = 1j/Re

 A[1][1-1] += a2 - ((4.0*a3)/yDelta**2)
 A[1][1]   += (a1*yDelta**2) - (2.0*a2) + ((6.0*a3)/yDelta**2)
 A[1][1+1] += a2 - ((4.0*a3)/yDelta**2)
 A[1][1+2] += a3/(yDelta**2)
 
 for i in range(2,numNodes-2):
  a1 = (-u[i]*alpha[k]**3) - (u2[i]*alpha[k]) + ((1j*alpha[k]**4)/Re)
  a2 = ( u[i]*alpha[k]) - ((2.0*1j*alpha[k]**2)/Re)
  a3 = 1j/Re
 
  A[i][i-2] += a3/(yDelta**2)
  A[i][i-1] += a2 - ((4.0*a3)/yDelta**2)
  A[i][i]   += (a1*yDelta**2) - (2.0*a2) + ((6.0*a3)/yDelta**2)
  A[i][i+1] += a2 - ((4.0*a3)/yDelta**2)
  A[i][i+2] += a3/(yDelta**2)
     

 a1 = (-u[numNodes-2]*alpha[k]**3) - (u2[numNodes-2]*alpha[k]) + ((1j*alpha[k]**4)/Re)
 a2 = ( u[numNodes-2]*alpha[k]) - ((2.0*1j*alpha[k]**2)/Re)
 a3 = 1j/Re

 A[numNodes-2][numNodes-2-2] += a3/(yDelta**2)
 A[numNodes-2][numNodes-2-1] += a2 - ((4.0*a3)/yDelta**2)
 A[numNodes-2][numNodes-2]   += (a1*yDelta**2) - (2.0*a2) + ((6.0*a3)/yDelta**2)
 A[numNodes-2][numNodes-2+1] += a2 - ((4.0*a3)/yDelta**2)
 
 A = np.delete(A,numNodes-1,0)
 A = np.delete(A,numNodes-1,1)
 A = np.delete(A,0,0)
 A = np.delete(A,0,1)
 

 B = np.zeros([numNodes,numNodes], dtype=complex)
 for i in range(1,numNodes-1):
  b1 = -alpha[k]**3
  b2 = alpha[k]
  B[i][i-1] += b2
  B[i][i]   += (b1*yDelta**2) - (2.0*b2)
  B[i][i+1] += b2

 B = np.delete(B,numNodes-1,0)
 B = np.delete(B,numNodes-1,1)
 B = np.delete(B,0,0)
 B = np.delete(B,0,1)
 


 # Pontes notes     
#  A[i][i-2] = (1.0/(yDelta**4))*(1.0/(alpha[k]*Re))*1j
#  A[i][i-1] = ((-4.0/(yDelta**4)) - (2.0*alpha[k]**2)*(1.0/(yDelta**2)))*(1.0/(alpha[k]*Re))*1j + u[i]*(1.0/(yDelta**2))
#  A[i][i]   = ((6.0/(yDelta**4)) - (2.0*alpha[k]**2)*(-2.0/(yDelta**2)) + (alpha[k]**4))*(1.0/(alpha[k]*Re))*1j + u[i]*(-2.0/(yDelta**2) - (alpha[k]**2)) - u2[i]
#  A[i][i+1] = ((-4.0/(yDelta**4)) - (2.0*alpha[k]**2)*(1.0/(yDelta**2)))*(1.0/(alpha[k]*Re))*1j + u[i]*(1.0/(yDelta**2))
#  A[i][i+2] = (1.0/(yDelta**4))*(1.0/(alpha[k]*Re))*1j
# 
#  B[i][i-1] = (1.0/(yDelta**2))
#  B[i][i]   = (-2.0/(yDelta**2) - (alpha[k]**2))
#  B[i][i+1] = (1.0/(yDelta**2))
#
# A = np.delete(A,-2,0)
# A = np.delete(A,-1,0)
# A = np.delete(A,1,0)
# A = np.delete(A,0,0)
# A = np.delete(A,-2,1)
# A = np.delete(A,-1,1)
# A = np.delete(A,1,1)
# A = np.delete(A,0,1)
# 
# B = np.delete(B,-2,0)
# B = np.delete(B,-1,0)
# B = np.delete(B,1,0)
# B = np.delete(B,0,0)
# B = np.delete(B,-2,1)
# B = np.delete(B,-1,1)
# B = np.delete(B,1,1)
# B = np.delete(B,0,1)


 invB = np.linalg.inv(B)
 eigenValue = np.linalg.eig(invB*A)
 c[k] = eigenValue[0].imag.max()


 #print alpha[k]
 print eigenValue
 #print ""
print c 
# --------------------------------------------------------------------------------------------
