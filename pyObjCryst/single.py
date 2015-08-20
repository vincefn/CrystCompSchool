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
plot(d.GetSinThetaOverLambda(),d.GetIobs()-d.GetIcalc())
