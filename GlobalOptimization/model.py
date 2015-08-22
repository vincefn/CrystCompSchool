import numpy as np
from scipy import weave
import time
import os

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

def Calc4(x,y,z,h,k,l):
  nhkl=len(h)
  natoms=len(x)
  fhklreal=np.zeros(nhkl)
  fhklimag=np.zeros(nhkl)
  err=weave.inline(  
        """
        for(int i=0;i<nhkl;i++)
          for(int i=0;i<natoms;i++)
          {
             const float tmp=6.28318530717958647692*(h(i)*x(i)+k(i)*y(i)+l(i)*z(i));
             fhklreal(i)+=cos(tmp);
             fhklimag(i)+=cos(tmp);
          }
        """,   
        ['fhklreal','fhklimag','x','y','z','h','k','l','natoms','nhkl'],   
        type_converters=weave.converters.blitz, 
         extra_compile_args=["-O3 -w -ffast-math -msse -msse2  -msse3 -msse4.1 -march=native -mfpmath=sse -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops -ftree-vectorize"],
        compiler="gcc"  
    )
  return fhklreal**2+fhklimag**2

  


def Calc5(x,y,z,h,k,l):
    # Using Weave and SSE - needs nbatoms to be a multiple of 4
    # requires sse_mathfun.h header in the directory
    code_CPU_fhkl_xyz="""
      const float PI2         = 6.28318530717958647692528676655900577f;
      for(unsigned long i=0;i<nhkl;i++)
      {
          float fr=0,fi=0;
          const float h=vh[i]*PI2;
          const float k=vk[i]*PI2;
          const float l=vl[i]*PI2;
          const float * __restrict__ px=vx;
          const float * __restrict__ py=vy;
          const float * __restrict__ pz=vz;
          __m128 vfr,vfi,vs,vc;
          float tmp[4];
          for(unsigned long at=0;at<natoms;at+=4)
          {
            float * __restrict__ ptmp=&tmp[0];
        
            // Dangerous ? Order of operation is not guaranteed - but it works...
            sincos_ps(_mm_set_ps(h* *px++ +k * *py++ + l * *pz++,
                                 h* *px++ +k * *py++ + l * *pz++,
                                 h* *px++ +k * *py++ + l * *pz++,
                                 h* *px++ +k * *py++ + l * *pz++),&vs,&vc);
            if(at==0) 
            {vfr=vc;vfi=vs;}
            else 
            {vfr=_mm_add_ps(vfr,vc);vfi=_mm_add_ps(vfi,vs);}
          }
          float tmp2[4];
          _mm_store_ps(tmp2,vfr);
          for(unsigned int j=0;j<4;++j) fr+=tmp2[j];
          _mm_store_ps(tmp2,vfi);
          for(unsigned int j=0;j<4;++j) fi+=tmp2[j];
          freal[i]=fr;
          fimag[i]=fi;
      }
    """
    nhkl=len(h)
    natoms=len(x)
    vh,vk,vl=h.astype(np.float32),k.astype(np.float32),l.astype(np.float32)
    vx,vy,vz=x.astype(np.float32),y.astype(np.float32),z.astype(np.float32)
    freal=np.zeros(nhkl,dtype=np.float32)
    fimag=np.zeros(nhkl,dtype=np.float32)
    pth= os.path.dirname(os.path.abspath(__file__))
    fhkl_xyz = weave.inline("Py_BEGIN_ALLOW_THREADS\n" + code_CPU_fhkl_xyz+ "Py_END_ALLOW_THREADS\n",
                                ['freal','fimag','vh', 'vk', 'vl', 'vx', 'vy', 'vz', 'nhkl', 'natoms'],
                                extra_compile_args=["-O3 -I"+pth+" -w -ffast-math -msse -msse2  -msse3 -msse4.1 -march=native -mfpmath=sse -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops -ftree-vectorize"],
                                compiler = 'gcc',
                                support_code="""#define USE_SSE2
                                                #include "sse_mathfun.h"
                                              """,
                                include_dirs=['./'])
    return freal**2+fimag**2

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

  n4=100
  dt4=time.time()
  for i in xrange(n4):
    iobs=Calc4(x,y,z,h,k,l)
  dt4=time.time()-dt4

  n5=100
  if (os.listdir('./').count('sse_mathfun.h')>0) and ((nbatoms%4)==0):
      dt5=time.time()
      for i in xrange(n5):
        iobs=Calc5(x,y,z,h,k,l)
      dt5=time.time()-dt5
  else:
      dt5=1e6

  print("dt1= %7.4fs  dt2= %7.4fs  dt3= %7.4fs  dt4= %7.4fs  dt5= %7.4fs"%(dt1,dt2,dt3,dt4,dt5))
  print("Calc1: %5.1f calc/s  Calc2: %5.1f calc/s  Calc3: %5.1f calc/s  Calc4: %5.1f calc/s  Calc5: %5.1f calc/s"%(n1/dt1,n2/dt2,n3/dt3,n4/dt4,n5/dt5))
  print("Calc1: %8.2f krefl.atoms/s  Calc2: %8.2f krefl.atoms/s  Calc3: %8.2f krefl.atoms/s  Calc4: %8.2f krefl.atoms/s  Calc5: %8.2f krefl.atoms/s"%(n1*nhkl/1e3/dt1,n2*nhkl/1e3/dt2,n3*nhkl/1e3/dt3,n4*nhkl/1e3/dt4,n5*nhkl/1e3/dt5))

testmodel(nbatoms=8,nh=3)
print()
testmodel(nbatoms=36,nh=10)
