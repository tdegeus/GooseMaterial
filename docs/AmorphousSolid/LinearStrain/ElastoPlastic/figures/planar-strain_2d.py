import numpy as np
import matplotlib.pyplot as plt
import GooseTensor as gt

plt.style.use(['goose','goose-latex'])

# rotation matrix, based on angle
rotation = lambda theta : np.array([
    [  np.cos(theta), -np.sin(theta) ],
    [  np.sin(theta),  np.cos(theta) ],
  ])

# simple shear strain
simple = np.array([
  [ 0. , 1. ],
  [ 1. , 0. ],
])

# normal
normal = np.array([[0.,1.]]).T

theta  = np.linspace(-np.pi/2.,+np.pi/2,800)
eps_ps = np.empty(theta.shape)
eps_ss = np.empty(theta.shape)
eps_eq = np.empty(theta.shape)

for i in range(theta.size):

  # apply a certain rotation
  R     = rotation(theta[i])
  eps   = (R.dot(simple)).dot(R.T)
  # n     = R.dot(normal)
  n = normal

  # get the planar strain
  s     = eps.dot(n)
  s    /= np.linalg.norm(s)
  sp    = s - s.T.dot(n) * n
  sp   /= np.linalg.norm(sp)
  gamma = (sp.T.dot(eps)).dot(n)
  epspn = gamma * ( sp.T*n + (sp.T*n).T )
  epss  = eps - epspn

  eps_eq[i] = np.sqrt(.5*gt.ddot22(eps  ,eps  ))
  eps_ss[i] = np.sqrt(.5*gt.ddot22(epspn,epspn))
  eps_ps[i] = np.sqrt(.5*gt.ddot22(epss ,epss ))

fig,ax = plt.subplots()

ax.plot(theta,eps_eq,label=r'$\varepsilon_\mathrm{eq}$')
ax.plot(theta,eps_ss,label=r'$\varepsilon_\mathrm{ss}$')
ax.plot(theta,eps_ps,label=r'$\varepsilon_\mathrm{ps}$')

plt.legend()

ax.xaxis.set_ticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'])
ax.xaxis.set_ticks([-np.pi/2.,-np.pi/4.,0,np.pi/4.,np.pi/2.])
ax.yaxis.set_ticks([-1,0,1])

plt.xlabel(r'$\theta$')
plt.ylabel(r'$\varepsilon$')

plt.show()
