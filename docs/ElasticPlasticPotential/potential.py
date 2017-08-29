
import numpy             as np
import matplotlib        as mpl
import matplotlib.pyplot as plt
plt.style.use(['goose','goose-latex','goose-tick-lower'])

# --------------------------------------------------------------------------------------------------

eps_m = np.linspace(-1.1,+1.1,1000)
U     = 9./2. * eps_m**2.

fig,ax = plt.subplots()

ax.plot(
  eps_m ,
  U     ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{m}$')

plt.ylabel(r'$U ( \varepsilon_\mathrm{m} ) $')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([])

plt.savefig('potential_U.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
V      = 3./2. * eps_eq**2.

fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  V      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$V ( \varepsilon_\mathrm{eq} ) $')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([])

plt.savefig('potential_V-elas.svg')


# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
V      = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  Vi  = .5 * ( (eps_eq-abar)**2. - delta**2. )

  V[ idx ] = Vi[ idx ]


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  V      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$V ( \varepsilon_\mathrm{eq} ) $')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([])

plt.savefig('potential_V-plas_raw.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
V      = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  Vi  = - ( delta / ( np.pi ) )**2. * ( 1. + np.cos( np.pi * ( eps_eq - abar ) / delta ) )

  V[ idx ] = Vi[ idx ]


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  V      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$V ( \varepsilon_\mathrm{eq} ) $')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([])

plt.savefig('potential_V-plas-smooth_raw.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
sig_eq = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  dVi = eps_eq-abar

  sig_eq[ idx ] = dVi[ idx ]


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  sig_eq      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\partial V / \partial \varepsilon_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('potential_dV-plas.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
sig_eq = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  dVi = ( delta / ( np.pi ) ) * np.sin( np.pi * ( eps_eq - abar ) / delta )

  sig_eq[ idx ] = dVi[ idx ]


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  sig_eq      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\partial V / \partial \varepsilon_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('potential_dV-plas-smooth.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
sig_eq = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  dVi = eps_eq-abar

  sig_eq[ idx ] = np.abs(dVi[ idx ])


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  sig_eq      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\sigma_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('potential_sigeq-plas.svg')

# --------------------------------------------------------------------------------------------------

eps_eq = np.linspace(0,10,1000)
sig_eq = np.zeros(eps_eq.shape)

a = [ -1. , 1. , 1.5 , 3. , 6. , 10.1 ]

for a0 , a1 in zip( a[:-1] , a[1:] ):

  abar  = ( a1 + a0 ) / 2.
  delta = ( a1 - a0 ) / 2.

  idx = np.where( ( eps_eq >= a0 ) * ( eps_eq  < a1 ) )[0]
  dVi = ( delta / ( np.pi ) ) * np.sin( np.pi * ( eps_eq - abar ) / delta )

  sig_eq[ idx ] = np.abs(dVi[ idx ])


fig,ax = plt.subplots()

ax.plot(
  eps_eq ,
  sig_eq      ,
  color = 'k'
)

plt.xlabel(r'$\varepsilon_\mathrm{eq}$')

plt.ylabel(r'$\sigma_\mathrm{eq}$')

ax.xaxis.set_ticks([0])
ax.yaxis.set_ticks([0])

plt.savefig('potential_sigeq-plas-smooth.svg')

# --------------------------------------------------------------------------------------------------

