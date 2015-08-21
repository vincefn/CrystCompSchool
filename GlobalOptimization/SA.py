from scipy import optimize
import model

nbatoms=10
x,y,z=model.genatoms(nbatoms)
h,k,l=model.genrefl(3)
nhkl=len(h)
