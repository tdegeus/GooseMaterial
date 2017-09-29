'''
Illustration of randomized behavior.

(c) T.W.J. de Geus | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid
'''

import numpy as np
import material_elastic_viscous as mat

import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex'])

# ==================================================================================================

eps  = 10.0
T    = 100.
ninc = 100000
dt   = T/float(ninc)

epsdot = (eps/T)*np.array([
  [0. , 1. , 0.],
  [1. , 0. , 0.],
  [0. , 0. , 0.],
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

# generate an exponential distribution of yield stresses: one for each increment
# (the model only select it at the instance it changes from elastic to plastic)
sigy_mean = 0.1
sigy_low  = 0.01
lam       = 1./(sigy_mean-sigy_low)
U         = np.random.random(ninc)
sigy      = sigy_low-np.log(1.-U)/lam

for inc in range(ninc):

  model.setNext_sigy(sigy[inc])
  sigma = model.stress(epsdot,dt,float(inc)*dt)
  model.increment()

  tau  [inc] = sigma [0,1]
  gamma[inc] = epsdot[0,1]*float(inc)*dt

fig,ax = plt.subplots()

ax.plot(gamma,tau*np.sqrt(3.))

ax.xaxis.set_ticks(np.linspace(0,10,6))
ax.yaxis.set_ticks(np.linspace(0,.3,4))

plt.xlim([-.05,10.05])
plt.ylim([-.01,0.32])

plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\tau \sqrt{3}$')

plt.savefig('example_random.svg')
plt.show()
