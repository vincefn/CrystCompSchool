from scipy import optimize
import numpy as np
import model
from pylab import plot,legend

# Warning: this only works with scipy before 0.16 (1 year old): optimize.anneal has now been removed from scipy !

nbatoms=10
x,y,z=model.genatoms(nbatoms)
h,k,l=model.genrefl(3)
nhkl=len(h)

scale=10
iobs=model.Calc3(x,y,z,h,k,l)*scale
sig=np.sqrt(iobs)
weight=1/iobs
iobs0=iobs.copy()
#iobs+=np.random.normal(0,1,iobs.shape)*iobs
iobs=np.random.poisson(iobs)

def chi2(xyz,doplot=False):
  n=len(xyz)//3
  icalc=model.Calc3(xyz[0:n],xyz[n:n*2],xyz[n*2:n*3],h,k,l)*scale
  c=(weight*(iobs-icalc)**2).sum()
  print(c)
  if doplot:
    clf()
    plot(iobs)
    plot(icalc)
    legen(('Obs','Calc'))
  return c

xyz=np.append(np.append(x,y),z)
n=len(xyz)//3
icalc=scale*model.Calc3(xyz[0:n],xyz[n:n*2],xyz[n*2:n*3],h,k,l)

xyz0=np.random.uniform(0,1,len(x)*3)

#xyz1=optimize.anneal(chi2,xyz0)
res=optimize.minimize(chi2,xyz0,method='Nelder-Mead')


