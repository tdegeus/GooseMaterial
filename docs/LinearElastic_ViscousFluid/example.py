
import numpy as np
import material_elastic_viscous as mat

import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex'])

fig = plt.figure(figsize=(21,6))
fig.set_tight_layout(True)

# ==================================================================================================

eps  = 1.0
T    = 10.
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1.],
  [1. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = mat.model(
  K      = 1.0,
  mu     = 1.0,
  sigy   = 0.1,
  tdamp  = 0.2,
  tfluid = 2.0,
  nd     = 2
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt,float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,1)

ax.plot(gamma,tau)

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau$')

# ==================================================================================================

eps  = 1.0
T    = 2.0
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1.],
  [1. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = mat.model(
  K      = 1.0,
  mu     = 1.0,
  sigy   = 0.1,
  tdamp  = 0.2,
  tfluid = 2.0,
  nd     = 2
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt,float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,2)

ax.plot(gamma,tau)

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau$')

# ==================================================================================================

eps  = 1.0
T    = 1.8
ninc = 40000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1.],
  [1. , 0.],
])

tau   = np.zeros((ninc))
gamma = np.zeros((ninc))

model = mat.model(
  K      = 1.0,
  mu     = 1.0,
  sigy   = 0.1,
  tdamp  = 0.2,
  tfluid = 2.0,
  nd     = 2
)

for inc in range(ninc):

  sigma = model.stress(epsdot,dt,float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

ax = fig.add_subplot(1,3,3)

ax.plot(gamma,tau)

ax.xaxis.set_ticks(np.linspace(0, 1,6))
ax.yaxis.set_ticks(np.linspace(0,.1,6))

plt.xlim([-.05,1.05])
plt.ylim([-.01,0.12])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau$')

plt.savefig('example.svg')
plt.show()

