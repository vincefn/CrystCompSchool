import numpy as np
import time

def genatoms(nbatoms):
# Random atomic positions (fractional coordinates)
  return np.random.uniform(0,1,(nbatoms)),np.random.uniform(0,1,(nbatoms)),np.random.uniform(0,1,(nbatoms))


def genrefl(nh):
  # Reflections
  h=np.arange(nh)
  k=np.arange(nh)[:,np.newaxis]
  l=np.arange(nh)[:,np.newaxis,np.newaxis]
  # Flatten the hkl arrays, remove (0,0,0)
  if True:
    h,k,l=(h+(k+l)*0).flatten()[1:], (k+(h+l)*0).flatten()[1:], (l+(h+k)*0).flatten()[1:]
  else:
    h=h+(k+l)*0
    k=k+(h+l)*0
    l=l+(h+k)*0
    h=h.flatten()[1:]
    k=k.flatten()[1:]
    l=l.flatten()[1:]
  return h,k,l

# Function to compute diffracted intensities
def Calc1(x,y,z,h,k,l):
  #Bad - two loops
  fhkl=np.zeros(h.shape,np.complex)
  for i in xrange(len(h)):
    for j in xrange(len(x)):
      fhkl[i]+=np.exp(2j*np.pi*(x[j]*h[i]+y[j]*k[i]+z[j]*l[i]))
  return abs(fhkl)**2

def Calc2(x,y,z,h,k,l):
  #Better - one vectorized loop
  fhkl=np.zeros(h.shape,np.complex)
  for i in xrange(len(h)):
    fhkl[i]=np.exp(2j*np.pi*(x*h[i]+y*k[i]+z*l[i])).sum()
  return abs(fhkl)**2

def Calc3(x,y,z,h,k,l):
  #Best (*) - no loop
  fhkl=np.exp(2j*np.pi*(x*h[:,np.newaxis]+y*k[:,np.newaxis]+z*l[:,np.newaxis])).sum(axis=1)
  return abs(fhkl)**2


# Try with 10 atoms, 3**3-1 reflections, and 50 atoms, 6**3-1 reflections
def testmodel(nbatoms=10,nh=3):
  x,y,z=genatoms(nbatoms)
  nbatoms=len(x)
  h,k,l=genrefl(nh)
  nhkl=len(h)

  
  n1=int(100/(nbatoms*nhkl/200))
  if n1<1:n1=1
  print n1
  dt1=time.time()
  for i in xrange(n1):
    iobs=Calc1(x,y,z,h,k,l)
  dt1=time.time()-dt1

  n2=100
  dt2=time.time()
  for i in xrange(n2):
    iobs=Calc2(x,y,z,h,k,l)
  dt2=time.time()-dt2

  n3=100
  dt3=time.time()
  for i in xrange(n3):
    iobs=Calc3(x,y,z,h,k,l)
  dt3=time.time()-dt3

  print("dt1= %7.4fs  dt1= %7.4fs  dt1= %7.4fs"%(dt1,dt2,dt3))
  print("Calc1: %5.1f calc/s  Calc2: %5.1f calc/s  Calc3: %5.1f calc/s"%(n1/dt1,n2/dt2,n3/dt3))
  print("Calc1: %8.2f krefl.atoms/s  Calc2: %8.2f krefl.atoms/s  Calc3: %8.2f krefl.atoms/s"%(n1*nhkl/1e3/dt1,n2*nhkl/1e3/dt2,n3*nhkl/1e3/dt3))


