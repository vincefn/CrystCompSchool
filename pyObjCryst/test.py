# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:25:35 2015

@author: vincent
"""

import pyobjcryst
import numpy as np
import platform
import os

#Get the data from IUCr journals (once)
if True:
  WGET="wget -q -O-"
  if platform.uname()[0]=='Darwin': WGET="curl "
  # http://scripts.iucr.org/cgi-bin/paper?dt3034
  os.system("%s 'http://scripts.iucr.org/cgi-bin/sendcif?dt3034sup1' > dt3034.cif"%WGET)
  os.system("%s 'http://journals.iucr.org/c/issues/2015/09/00/dt3034/dt3034Isup2.hkl' > dt3034Isup2.hkl"%WGET)
  # http://scripts.iucr.org/cgi-bin/paper?hx5075
  os.system("%s 'http://scripts.iucr.org/cgi-bin/sendcif?hx5075sup1' > hx5075.cif"%WGET)
  os.system("%s 'http://scripts.iucr.org/cgi-bin/sendsupfiles?hx5075&file=hx5075sup9.txt' > hx5075sup9.txt"%WGET)

c=pyobjcryst.crystal.CreateCrystalFromCIF(open("dt3034.cif"),True,True)

d=pyobjcryst.diffractiondatasinglecrystal.DiffractionDataSingleCrystal(c)
d.ImportCIF("dt3034Isup2.hkl")
d.PrintObsCalcData()
print(d.GetRw())
d.SetWavelength(0.71073)
print(d.GetRw())

#print(c.GetNbScatterer())
#c.ConnectAtoms()
print(c.GetNbScatterer())
m=c.GetScatterer(0)
print(m.GetFormula())
#plot(d.GetSinThetaOverLambda(),d.GetIobs()-d.GetIcalc())


###################

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
from pylab import plot,pi
plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternObs())
plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternCalc())

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

#from pylab import plot,pi
plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternObs())
plot(p.GetPowderPatternX()*180/pi,p.GetPowderPatternCalc())

pd.SetExtractionMode(True,True)
pd.ExtractLeBail(10)

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

lsqr.FixAllPar()
lsqr.Print()
lsq.SetParIsFixed("Zero",False)
lsq.SetParIsFixed("2ThetaDispl",False)
lsq.SetParIsFixed("2ThetaTransp",False)
lsq.Refine(20,True)

lsq.SetParIsFixed("a",False)
lsq.SetParIsFixed("b",False)
lsq.SetParIsFixed("c",False)
lsq.SetParIsFixed("alpha",False)
lsq.SetParIsFixed("beta",False)
lsq.SetParIsFixed("gamma",False)
lsq.SetParIsFixed("W",False)
lsq.Refine(20,True)
lsq.Refine(20,True)
lsq.SetParIsFixed("U",False)
lsq.SetParIsFixed("V",False)
lsq.Refine(20,True)
lsq.Refine(20,True)
lsq.SetParIsFixed("Eta0",False)
lsq.Refine(20,True)
lsq.SetParIsFixed("Eta1",False)
lsq.Refine(20,True)
lsq.SetParIsFixed("Asym0",False)
lsq.SetParIsFixed("Asym1",False)
lsq.Refine(20,True)

b.GetPar(0)
b.GetPar(0).GetType()
t=b.GetPar(0).GetType()
t.GetName()
lsq.SetParIsFixed(t,False)
b.FixParametersBeyondMaxresolution(lsq)
lsqr.Print()
lsq.Refine(20,True)

pd.ExtractLeBail(10)
lsq.Refine(20,True)


##############################
pd.SetExtractionMode(False)
c.RandomizeConfiguration()

m=pyobjcryst.montecarloobj.MonteCarloObj()
m.AddRefinableObj(c)
m.AddRefinableObj(p)
m.Optimize(100000)


