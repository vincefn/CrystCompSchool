import model
import numpy
from pylab import *

def PT(p,func_chi2,nbtrial=1e5,nbworld=10,cost_history={}):
  if len(cost_history.keys())>0:
    history_shift=max(cost_history[0][0])
  else:
    history_shift=0
    #For each world we keep 1) trial # 2) chi2
    cost_history['best']=[[],[]]
    for i in xrange(nbworld):cost_history[i]=[[],[]]
  nbtrial0=nbtrial
  nbpar=len(p)
  temperature=linspace(100,10000,nbworld).astype(float32)
  #dpar=0.02 # Base displacement amplitude for parameters
  amplitude_min=0.0002
  amplitude_max=0.25
  amplitude=amplitude_min*(amplitude_max/amplitude_min)**(arange(nbworld,dtype=float32)/(nbworld-1))
  nbaccept=zeros(nbworld)
  acceptance_rate=zeros(nbworld)
  bestcost=func_chi2(p)
  bestpar=p.copy()
  cost=ones(nbworld)*bestcost
  p_world={}
  for world in xrange(nbworld):
    if world==0:
      p_world[world]=p.copy()
    else:
      p_world[world]=numpy.random.uniform(0,1,len(p))
  p1=p.copy()
  nbTestPerWorld=200
  while nbtrial>0:
    for world in xrange(nbworld):
      nbaccept[world]=0
      for i in xrange(nbTestPerWorld):
        nbtrial-=1
        lastcost=cost[world]
        p1=p_world[world]+uniform(-1,1,p1.shape)*amplitude[world]
        #p1=p_world[world]+uniform(-1,1,p1.shape)*amplitude[world]*(uniform(0,1,p1.shape)>0.7)
        newcost=func_chi2(p1)
        if newcost<lastcost:
          acc=1
          p_world[world]=p1.copy()
          lastcost=newcost
          cost[world]=newcost
          if newcost<bestcost:
            bestcost=newcost
            bestpar=p1.copy()
            #print("Trial #%8d-NEW best cost, world #%2d, Chi^2=%12.2f"%(nbtrial,world,newcost))
            cost_history['best'][0].append(nbtrial0-nbtrial+history_shift)
            cost_history['best'][1].append(bestcost)
          nbaccept[world]+=1
        elif exp((lastcost-newcost)/temperature[world])>uniform(0,1):
          cost[world]=newcost
          p_world[world]=p1.copy()
          nbaccept[world]+=1
          acc=1
        else:acc=0
        #print("World #%2d, current cost=%12.2f, accept=%d  (best cost=%12.2f)"%(world,newcost,acc,bestcost))
      acceptance_rate[world]=nbaccept[world]/float(nbTestPerWorld)
      #print("World #%2d, acceptance rate=%3.0f%%, current cost=%12.2f, T=%8.2f, amplitude=%8.5f"%(world,acceptance_rate[world]*100,cost[world],temperature[world],amplitude[world]))
      if acceptance_rate[world]<.1:
        temperature[world]*=1.5
      elif acceptance_rate[world]>.3:
        temperature[world]/=1.5
    for world in xrange(nbworld-1,0,-1):
      if exp((cost[world-1]-cost[world])/temperature[world-1])>uniform(0,1):
        p_world[world],p_world[world-1]=(p_world[world-1]).copy(),(p_world[world]).copy()
        cost[world]   ,cost[world-1]   =cost[world-1]             ,cost[world]
        print("Exchange: %d <-> %d"%(world,world-1))
    for world in xrange(nbworld-1,-1,-1):
      print("World #%2d, acceptance rate=%3.0f%%, current cost=%12.2f, T=%12.8f, amplitude=%8.5f"%(world,acceptance_rate[world]*100,cost[world],temperature[world],amplitude[world]))
      cost_history[world][0].append(nbtrial0-nbtrial+history_shift)
      cost_history[world][1].append(cost[world])
    print("Trial #%8d, current best cost=%12.2f"%(nbtrial,bestcost))
    
  return bestpar



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
xyz1=PT(xyz0,chi2,nbtrial=3e5,nbworld=10,cost_history=cost_history)

x,y=cost_history['best']
semilogy(x,y)
xyz1=PT(xyz1,chi2,nbtrial=5e5,nbworld=30,cost_history=cost_history)
