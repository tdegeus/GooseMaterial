'''
Material response to simple-shear deformation. Three different loading rates are included to show
that the 'fluid' can soften or harden depending on the ratio of the loading and damping rate
compared to the ratio of the yield stress and the shear modulus.

See "notes_analytical-simple-shear.pdf"

(c) T.W.J. de Geus | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid
'''

import numpy as np
import GooseSolid as gs

import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex'])

fig = plt.figure(figsize=(21,6))
fig.set_tight_layout(True)

# ==================================================================================================

K      = 1.0
G      = 1.0
sigy   = 0.1
Tdamp  = 0.1
Tfluid = 2.0

# ==================================================================================================

eps  = 1.0
T    = 5./( sigy/np.sqrt(3.) / Tdamp / (2.*G) )
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = gs.LinearElastic_ViscousFluid(
  K      = K     ,
  G      = G     ,
  sigy   = sigy  ,
  Tdamp  = Tdamp ,
  Tfluid = Tfluid
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,1)

ax.plot(gamma,tau*np.sqrt(3.))

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')

plt.title(r'$T_e = 5 \frac{2 G}{\sigma_\mathrm{y}/\sqrt{3}} \; T$')

# ==================================================================================================

eps  = 1.0
T    = 1./( sigy/np.sqrt(3.) / Tdamp / (2.*G) )
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = gs.LinearElastic_ViscousFluid(
  K      = K     ,
  G      = G     ,
  sigy   = sigy  ,
  Tdamp  = Tdamp ,
  Tfluid = Tfluid
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,2)

ax.plot(gamma,tau*np.sqrt(3.))

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')

plt.title(r'$T_e = \frac{2 G}{\sigma_\mathrm{y}/\sqrt{3}} \; T$')

# ==================================================================================================

eps  = 1.0
T    = .9/( sigy/np.sqrt(3.) / Tdamp / (2.*G) )
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = gs.LinearElastic_ViscousFluid(
  K      = K     ,
  G      = G     ,
  sigy   = sigy  ,
  Tdamp  = Tdamp ,
  Tfluid = Tfluid
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,3)

ax.plot(gamma,tau*np.sqrt(3.))

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')

plt.title(r'$T_e = \frac{9}{10} \frac{2 G}{\sigma_\mathrm{y}/\sqrt{3}} \; T$')

plt.savefig('example.svg')
plt.show()

