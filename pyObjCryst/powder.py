# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:25:35 2015

@author: vincent
"""

import pyobjcryst
import numpy as np
import platform
import os
from pylab import plot,pi,figure,legend,xlabel,ylabel,xlim,ylim,clf

def pshow(p):
    figure(1)
    # Try to keep previous limits
    x0,x1=xlim()
    y0,y1=ylim()
    clf()
    plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternObs(),'b')
    plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternCalc(),'r')
    legend(('Obs','Calc'))
    xlabel('2Theta')
    ylabel('Intensity')
    if x1>2:
        xlim(x0,x1)
        ylim(y0,y1)
    draw()

#Get the data from IUCr journals (once)
if False:
  WGET="wget -q -O-"
  if platform.uname()[0]=='Darwin': WGET="curl "
  # http://scripts.iucr.org/cgi-bin/paper?hx5075
  os.system("%s 'http://scripts.iucr.org/cgi-bin/sendcif?hx5075sup1' > hx5075.cif"%WGET)
  os.system("%s 'http://scripts.iucr.org/cgi-bin/sendsupfiles?hx5075&file=hx5075sup9.txt' > hx5075sup9.txt"%WGET)

c=pyobjcryst.crystal.CreateCrystalFromCIF(open("hx5075.cif"),True,True)
c.GetScatterer(0).GetFormula()
p=pyobjcryst.powderpattern.PowderPattern()
p.ImportPowderPattern2ThetaObsSigma("hx5075sup9.txt")
x=p.GetPowderPatternX()

# Background
bx=np.linspace(x.min(),x.max(),20)
by=np.zeros(bx.shape)

b=pyobjcryst.powderpatternbackground.PowderPatternBackground()
b.SetInterpPoints(bx,by)
p.AddPowderPatternComponent(b)
b.Print()
b.UnFixAllPar()
b.OptimizeBayesianBackground()

pshow(p)

pd=pyobjcryst.powderpatterndiffraction.PowderPatternDiffraction()
pd.SetCrystal(c)

c.GetNbOption()
c.GetOption(0).GetName()
c.GetOption(1).GetName()
c.GetOption(1).GetChoiceName(0)
c.GetOption(1).GetChoiceName(1)
c.GetOption(1).SetChoice(0)


p.AddPowderPatternComponent(pd)
p.SetWavelength(0.801203)
p.Prepare()
pd.SetReflectionProfilePar(pyobjcryst.powderpatterndiffraction.ReflectionProfileType.PROFILE_PSEUDO_VOIGT,0.0000001)
p.FitScaleFactorForIntegratedRw()

pshow(p)


# LSQ
lsq=pyobjcryst.lsqnumobj.LSQNumObj()
lsq.SetRefinedObj(p)
lsq.GetCompiledRefinedObj().Print()
lsq.PrepareRefParList(True)
lsq.GetCompiledRefinedObj().Print()
lsq.SetRefinedObj(p,0,True,True)
lsq.PrepareRefParList(True)
lsq.GetCompiledRefinedObj().Print()
lsqr=lsq.GetCompiledRefinedObj()

p.SetMaxSinThetaOvLambda(0.5)

# Le Bail
pd.SetExtractionMode(True,True)
pd.ExtractLeBail(20)
pshow(p)

for i in range(4):
    # Profile fitting
    lsqr.FixAllPar()
    lsqr.Print()
    lsq.SetParIsFixed("Zero",False)
    lsq.SetParIsFixed("W",False)
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)

    lsq.SetParIsFixed("U",False)
    lsq.SetParIsFixed("V",False)
    lsq.SetParIsFixed("Eta0",False)
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)

    lsq.SetParIsFixed("Eta1",False)
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)

    lsq.SetParIsFixed("Asym0",False)
    lsq.SetParIsFixed("Asym1",False)
    lsq.SetParIsFixed("Asym2",False)
    lsq.SetParIsFixed("2ThetaDispl",False)
    lsq.SetParIsFixed("2ThetaTransp",False)
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)

    # Background
    b.GetPar(0)
    b.GetPar(0).GetType()
    t=b.GetPar(0).GetType()
    t.GetName()
    lsq.SetParIsFixed(t,False) # Unfix a type of parameters
    b.FixParametersBeyondMaxresolution(lsqr)
    lsqr.Print()
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)

    # Unit Cell
    lsq.SetParIsFixed("a",False)
    lsq.SetParIsFixed("b",False)
    lsq.SetParIsFixed("c",False)
    lsq.SetParIsFixed("alpha",False)
    lsq.SetParIsFixed("beta",False)
    lsq.SetParIsFixed("gamma",False)
    lsq.Refine(5,True)
    pd.ExtractLeBail(2)
    pshow(p)



###################################
pd.SetExtractionMode(False)
c.RandomizeConfiguration() # Don't cheat !

m=pyobjcryst.montecarloobj.MonteCarloObj()
m.AddRefinableObj(c)
m.AddRefinableObj(p)
m.Optimize(400000)
