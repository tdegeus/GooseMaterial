
import GooseTensor as gt
import GooseSolid  as gs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['goose','goose-latex','goose-tick-lower'])

# --------------------------------------------------------------------------------------------------

mat = gs.ElasticPlasticPotential( 1. , 1. , [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ] , False )

Eps = np.array([
  [ 0. , 1. , 0. ],
  [ 1. , 0. , 0. ],
  [ 0. , 0. , 0. ],
])

ninc   = 20000
eps_eq = np.zeros((ninc))
sig_eq = np.zeros((ninc))

for i,d in enumerate(np.linspace(0,10.*np.sqrt(3.)/2.,ninc)):

  eps       = d * Eps
  sig       = mat.stress(eps)
  sig_eq[i] = np.sqrt(3./2.*gt.ddot22(sig,sig))
  eps_eq[i] = np.sqrt(2./3.*gt.ddot22(eps,eps))

plt.style.use(['goose','goose-latex'])

fig,ax = plt.subplots()

ax.plot(eps_eq,sig_eq)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\sigma_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('stress-strain.svg')
plt.show()


# --------------------------------------------------------------------------------------------------

mat = gs.ElasticPlasticPotential( 1. , 1. , [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ] , True )

Eps = np.array([
  [ 0. , 1. , 0. ],
  [ 1. , 0. , 0. ],
  [ 0. , 0. , 0. ],
])

ninc   = 20000
eps_eq = np.zeros((ninc))
sig_eq = np.zeros((ninc))

for i,d in enumerate(np.linspace(0,10.*np.sqrt(3.)/2.,ninc)):

  eps       = d * Eps
  sig       = mat.stress(eps)
  sig_eq[i] = np.sqrt(3./2.*gt.ddot22(sig,sig))
  eps_eq[i] = np.sqrt(2./3.*gt.ddot22(eps,eps))

fig,ax = plt.subplots()

ax.plot(eps_eq,sig_eq)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\sigma_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('stress-strain-smooth.svg')
plt.show()


# --------------------------------------------------------------------------------------------------

