'''
Effect of different parameter variations for a pure-shear strain.

(c) T.W.J. de Geus | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid
'''

import numpy       as np
import GooseSolid  as gs
import GooseTensor as gt

import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex'])

# ==================================================================================================
# calculate stress strain response
# ==================================================================================================

def response(gamma0,n,sigy0,H,m,dt):

  # convert elastic constants
  K,G   = gs.ConvertElasticParameters("E,nu",1.0,0.3,"K,G")

  # define the material
  model = gs.ViscoPlasticHardeningLinearElastic(
    K      = K,
    G      = G,
    gamma0 = gamma0,
    n      = n     ,
    sigy0  = sigy0 ,
    H      = H     ,
    m      = m     ,
  )

  # define stain
  E    = 0.1 * np.sqrt(3.)/2.
  ninc = 100
  eps  = np.zeros((3,3))
  deps = np.array([
    [ 1. ,  0. ,  0.],
    [ 0. , -1. ,  0.],
    [ 0. ,  0. ,  0.],
  ])*E/float(ninc)

  eps_eq = np.zeros((ninc+1),dtype='float64')
  sig_eq = np.zeros((ninc+1),dtype='float64')

  for inc in range(ninc):

    eps += deps
    sig  = model.stress(eps,dt)

    model.increment()

    sig_d = sig - gt.trace2( sig )/3.

    sig_eq[ inc + 1 ] = np.sqrt( 3./2. * gt.ddot22( sig_d , sig_d ) )
    eps_eq[ inc + 1 ] = np.sqrt( 2./3. * gt.ddot22( eps   , eps   ) )

  return ( eps_eq , sig_eq )

# ==================================================================================================
# plot
# ==================================================================================================

# open figure
fig = plt.figure(figsize=(23,14))
fig.set_tight_layout(True)

# reference slip rate
# -------------------

ax   = fig.add_subplot(2,3,1)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([0.01,0.10,1.00,10.0,100.0]):

  parameters['gamma0'] = var
  eps_eq, sig_eq       = response(**parameters)
  col                  = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$\dot{\gamma}_0 = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

# hardening exponent
# ------------------

ax   = fig.add_subplot(2,3,2)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([0.05,0.1,0.2,0.4,0.8]):

  parameters['n'] = var
  eps_eq, sig_eq  = response(**parameters)
  col             = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$n = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

# slip rate
# ---------

ax   = fig.add_subplot(2,3,3)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([1.0e-5,1.0e-4,1.0e-3,1.0e-2,1.0e-1]):

  parameters['dt'] = var
  eps_eq, sig_eq   = response(**parameters)
  col              = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$\Delta t = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

# initial yield stress
# --------------------

ax   = fig.add_subplot(2,3,4)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([0.008,0.009,0.01,0.011,0.012]):

  parameters['sigy0'] = var
  eps_eq, sig_eq      = response(**parameters)
  col                 = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$\sigma_\mathrm{y0} = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

# hardening modulus
# -----------------

ax   = fig.add_subplot(2,3,5)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([0.005,0.01,0.02,0.04,0.08]):

  parameters['H'] = var
  eps_eq, sig_eq  = response(**parameters)
  col             = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$H = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

# hardening exponent
# ------------------

ax   = fig.add_subplot(2,3,6)

# set colormap
cmap = plt.get_cmap('RdBu_r',5)

# set parameters
parameters = dict(
  gamma0 = 1.00,
  n      = 0.20,
  sigy0  = 0.01,
  H      = 0.02,
  m      = 1.00,
  dt     = 1.0e-3,
)

# vary rate-exponent
for i,var in enumerate([0.6,0.8,1.0,1.2,1.4]):

  parameters['n'] = var
  eps_eq, sig_eq  = response(**parameters)
  col             = cmap(i) if i != 2 else 'k'

  ax.plot(eps_eq, sig_eq,linewidth=2.,color=col,label=r'$n = %s$'%str(var))

ax.xaxis.set_ticks([.00,.05,.10])
ax.yaxis.set_ticks([.00,.01,.02,.03])

plt.xlim([-.005,.105])
plt.ylim([-.001,.031])

plt.xlabel(r'$\bar{\varepsilon}$')
plt.ylabel(r'$\bar{\sigma}_\mathrm{eq} / E$')

plt.legend(loc='lower right')

plt.savefig('parameters.svg')
plt.show()
