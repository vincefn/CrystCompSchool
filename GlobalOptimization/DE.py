import model
import numpy
import time
from pylab import *
from scipy import optimize

def DE(nbgeneration,p0,func_chi2,cr=0.5,f=0.8,nbworld=500,strategy="DE/best/1/exp",prefix=None):
   nbpar=len(p0)
   nbgeneration0=nbgeneration
   worlds=[]
   worldstrial=[]
   bestworld=0
   rhomin=0.8
   rhomax=1.2
   newbest=False
   cost=zeros(nbworld)
   for i in xrange(nbworld):
      worlds.append(p0.copy())
      if i!=0: # Randomize all but one world
        worlds[i]=numpy.random.uniform(0,1,p0.shape)
      cost[i]=func_chi2(worlds[i])
      worldstrial.append(worlds[i].copy())
      print "DE: init world %3d: Chi2=%12.2f"%(i,cost[i])
   bestcost=cost.min()
   bestworld=argmin(cost)
   bestpar=worlds[bestworld]
   t0=time.time()
   while nbgeneration>0:
      for i in xrange(nbworld):
        done=False
        if strategy=="DE/best/1/exp":
          r1,r2=0,0
          while (r1==i) or (r1==bestworld):
              r1=numpy.random.randint(nbworld)
          while (r2==i) or (r2==bestworld) or(r2==r1):
              r2=numpy.random.randint(nbworld)
          crtmp=numpy.random.uniform(0,1,nbpar)<cr
          worldstrial[i]=worlds[i]*(crtmp==False)+(worlds[bestworld]+f*(worlds[r1]-worlds[r2]))*(crtmp==True)
          #print worldstrial[i]-worlds[i]
          done=True
        if strategy=="DE/rand/1/exp":
          r1,r2,r3=0,0,0
          while (r1==i) :
            r1=numpy.random.randint(nbworld)
          while (r2==i) or(r2==r1):
            r2=numpy.random.randint(nbworld)
          while (r3==i) or(r3==r1) or (r3==r2):
            r3=numpy.random.randint(nbworld)
          crtmp=numpy.random.uniform(0,1,nbpar)<cr
          worldstrial[i]=worlds[i]*(crtmp==False)+(worlds[r3]+f*(worlds[r1]-worlds[r2]))*(crtmp==True)
          done=True
        if strategy=="DE/rand-to-best/1/exp":
          r1,r2=0,0
          while (r1==i) or (r1==bestworld):
              r1=numpy.random.randint(nbworld)
          while (r2==i) or (r2==bestworld) or(r2==r1):
              r2=numpy.random.randint(nbworld)
          crtmp=numpy.random.uniform(0,1,nbpar)<cr
          worldstrial[i]=worlds[i]*(crtmp==False)+(worlds[i]+f*(worlds[bestworld]-worlds[i]+worlds[r1]-worlds[r2]))*(crtmp==True)
          done=True
        assert(done)
      for i in xrange(nbworld):
        c=func_chi2(worldstrial[i])
        #print "DE(%4d): Chi2(%3d)=%12.2f->%12.2f"%(nbgeneration,i,cost[i],c)
        if c<cost[i]:
          worlds[i]=worldstrial[i].copy()
          cost[i]=c
          if c<bestcost:
            bestcost=c
            bestpar=worlds[i]
            bestworld=i
            newbest=True
            print "NEW BEST - DE(%4d): Chi2(%3d)=%12.2f"%(nbgeneration,i,c)
      ## test for stagnation ?
      nbgeneration-=1
      dt=(time.time()-t0)/(nbgeneration0-nbgeneration)
      print "DE: generations to go:%d (%5.2fh remaining, %5.2fs/generation), best cost=%5.3f"%(nbgeneration,dt*nbgeneration/3600,dt,bestcost)
        
   return worlds[bestworld]

nbatoms=5
x,y,z=model.genatoms(nbatoms)
h,k,l=model.genrefl(5)
nhkl=len(h)

scale=10
iobs=model.Calc3(x,y,z,h,k,l)*scale
sig=sqrt(iobs)
weight=1/iobs
iobs0=iobs.copy()
#iobs+=np.random.normal(0,1,iobs.shape)*iobs
iobs=np.random.poisson(iobs)

def chi2(xyz,doplot=False):
  n=len(xyz)//3
  icalc=model.Calc3(xyz[0:n],xyz[n:n*2],xyz[n*2:n*3],h,k,l)*scale
  c=(weight*(iobs-icalc)**2).sum()
  #print(c)
  if doplot:
    clf()
    plot(iobs)
    plot(icalc)
    legen(('Obs','Calc'))
  return c


xyz0=np.random.uniform(0,1,len(x)*3)
cost_history={}
chi2(xyz0)
xyz1=DE(500,xyz0,chi2,cr=0.7,f=0.8,nbworld=500,strategy="DE/rand-to-best/1/exp")

#res=optimize.differential_evolution(chi2, [(0,1)]*len(xyz0), strategy='best1bin', maxiter=1000)
