'''
Example of linear and non-linear hardening.

(c) T.W.J. de Geus | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid
'''

import numpy as np
import GooseSolid as gs

import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex'])

fig = plt.figure(figsize=(14,6))
fig.set_tight_layout(True)

# ==================================================================================================

eps  = 1.0
T    = 1.
ninc = 1000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

K,G   = gs.ConvertElasticParameters("E,nu",1.0,0.3,"K,G")

model = gs.PlasticLinearElastic(
  K      = K,
  G      = G,
  sigy0  = 0.08,
  H      = 0.02,
)

for inc in range(ninc):

  sigma = model.stress(epsdot*float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,2,1)

ax.plot(gamma,tau*np.sqrt(3.0))

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')
plt.title (r'$m = 1$')

# ==================================================================================================

eps  = 1.0
T    = 1.
ninc = 1000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

K,G   = gs.ConvertElasticParameters("E,nu",1.0,0.3,"K,G")

model = gs.PlasticLinearElastic(
  K      = K,
  G      = G,
  sigy0  = 0.08,
  H      = 0.02,
  m      = 0.2 ,
)

for inc in range(ninc):

  sigma = model.stress(epsdot*float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,2,2)

ax.plot(gamma,tau*np.sqrt(3.0))

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')
plt.title (r'$m = 0.2$')

# ==================================================================================================

plt.savefig('example.svg')
plt.show()

