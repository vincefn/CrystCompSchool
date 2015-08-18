# CrystCompSchool
Tutorials for the crystallographic computing school


# Tutorial #1: (py)ObjCryst++ (1h30)
The ObjCryst++ C++ library (http://objcryst.sf.net) was developed for ab initio structure solution from powder diffraction. In this tutorial we will first give a quick tour of the C++ library capabilities. Then we will use the python interface to exploit this library, including: creating or importing crystal structure from CIF, calculating & displaying powder patterns, performing a global optimization.

Keywords: python/C++, small molecules/inorganic structures, powder diffraction
Requirements: install python+libobjcryst+pyobjcryst+matplotlib, optional install of Fox, good python knowledge


# Tutorial #2: global optimization (1h to 1h30 ?)
Finding the numerical solution to a problem can be tricky when the 'parameter space' to explore is large. In this tutorial we focus on 'well-defined' problems where the parametrization is known (e.g. a fixed number of atoms) but the value range of each parameter is very large. We will illustrate a few algorithms (simulated annealing, parallel tempering,differential evolution) which can be used to tackle this type of problems

Requirements: python/numpy/matplotlib, some python experience


# Tutorial #3: quick introduction to GPU computing with pyOpenCL and PyNX (40')
GPU use for general computing has been developed a lot during the past 10 years, due to the significant speedup (up to 1000x) it can provide vs traditional microprocessing computing. After a quick introduction about GPU computing, the notion of massively parallel calculations and CUDA/OpenCL, we will try simple examples using the PyNX library which allows fast calculations of structure factors.

Requirements: pyOpenCL (or pyCUDA), PyNX

# Tutorial #4: quick introduction to databases, SQL queries, and the Crystallographic Open Database (30')
Databases are often perceived as complicated by computing scientists. This quick tutorial will show that accessing (and even creating) databases is in fact trivial. We will give a very quick introduction to SQL queries, and how to search and retrieve data from the Crystallographic Open Database (www.crystallography.net/) as an example.


Requirements: python with mysql, sqlite modules, mysql, internet access. Beginners only