import numpy as np
from numpy import linalg

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


Re = 5772
alpha = np.zeros([8,1],dtype=float)
alpha[0] = 0.0
alpha[1] = 0.2
alpha[2] = 0.4
alpha[3] = 0.6
alpha[4] = 0.8
alpha[5] = 1.0
alpha[6] = 1.2
alpha[7] = 1.4
c = np.zeros([len(alpha),1],dtype=float)

for k in range(0,len(alpha)):
 A = np.zeros([numNodes,numNodes], dtype=complex)
 B = np.zeros([numNodes,numNodes], dtype=float)
 for i in range(2,numNodes-2):
  A[i][i-2] = (1.0/(yDelta**4))*(1.0/(alpha[k]*Re))*1j
  A[i][i-1] = ((-4.0/(yDelta**4)) - (2.0*alpha[k]**2)*(1.0/(yDelta**2)))*(1.0/(alpha[k]*Re))*1j + u[i]*(1.0/(yDelta**2))
  A[i][i]   = ((6.0/(yDelta**4)) - (2.0*alpha[k]**2)*(-2.0/(yDelta**2)) + (alpha[k]**4))*(1.0/(alpha[k]*Re))*1j + u[i]*(-2.0/(yDelta**2) - (alpha[k]**2)) - u2[i]
  A[i][i+1] = ((-4.0/(yDelta**4)) - (2.0*alpha[k]**2)*(1.0/(yDelta**2)))*(1.0/(alpha[k]*Re))*1j + u[i]*(1.0/(yDelta**2))
  A[i][i+2] = (1.0/(yDelta**4))*(1.0/(alpha[k]*Re))*1j
 
  B[i][i-1] = (1.0/(yDelta**2))
  B[i][i]   = (-2.0/(yDelta**2) - (alpha[k]**2))
  B[i][i+1] = (1.0/(yDelta**2))
 
 
 #A[0][0] = 1.0
 #A[0][1] = 1.0
 #A[1][0] = -1.0
 #A[1][1] = 1.0
 #A[numNodes-1][numNodes-1] = 1.0
 #A[numNodes-1][numNodes-2] = 1.0
 #A[numNodes-2][numNodes-1] = 1.0
 #A[numNodes-2][numNodes-2] = -1.0
 
 A = np.delete(A,-2,0)
 A = np.delete(A,-1,0)
 A = np.delete(A,1,0)
 A = np.delete(A,0,0)
 A = np.delete(A,-2,1)
 A = np.delete(A,-1,1)
 A = np.delete(A,1,1)
 A = np.delete(A,0,1)
 
 B = np.delete(B,-2,0)
 B = np.delete(B,-1,0)
 B = np.delete(B,1,0)
 B = np.delete(B,0,0)
 B = np.delete(B,-2,1)
 B = np.delete(B,-1,1)
 B = np.delete(B,1,1)
 B = np.delete(B,0,1)
 
 invB = np.linalg.inv(B)
 eigenValue = np.linalg.eig(np.dot(invB,A))
 c[k] = eigenValue[0].imag.max()
 
print c 
