
import numpy             as np
import GooseSolid        as gs
import GooseTensor       as gt
import matplotlib.pyplot as plt

plt.style.use(['goose','goose-latex'])

# ----------------------------------------- DEFINE MATERIAL ----------------------------------------

# define material model
mat = gs.NonLinearElastic(
  K     =  0.833333,
  sig0  =  0.5,
  eps0  =  0.1,
  n     = 10.1,
)

# ------------------------------------------- PRE-LOADING ------------------------------------------

# define strain tensor
eps       = np.zeros((3,3))
eps[0,1] += .1
eps[1,0] += .1

K4, sig = mat.tangent_stress( eps )

# ---------------------------------------- CONSISTENCY CHECK ---------------------------------------

# set perturbations, initialize the result
ndelta = 201
DELTA  = np.zeros((ndelta))
ETA    = np.zeros((ndelta))

# apply different perturbations, using a random strain
for i,delta in enumerate(np.logspace(-18,0,ndelta)):

  deps     = delta * ( np.random.random((3,3)) - 0.5 )
  deps     = .5 * ( deps + deps.T ) # NB has to be symmetrized because the module symmetrizes
  _,sig_n  = mat.tangent_stress( eps + deps )
  dsig     = sig_n - sig
  ETA  [i] = np.linalg.norm( dsig - gt.ddot42(K4,deps) ) / np.linalg.norm( dsig )
  DELTA[i] = np.linalg.norm( deps                      )

# ---------------------------------------------- PLOT ----------------------------------------------

# sort for intuitive plot
idx   = np.argsort(DELTA)
DELTA = DELTA[idx]
ETA   = ETA  [idx]

# plot the result
fig,ax = plt.subplots()

ax.plot([-18,0],[-18,  0],color='k',linestyle='--')
ax.plot([-18,0],[  0,-18],color='k',linestyle='--')
ax.plot(np.log(DELTA)/np.log(10.),np.log(ETA)/np.log(10.),color='r')

ax.set_title('Consistency check')

plt.xlim([-18,0])
plt.ylim([-18,0])

plt.xlabel(r'$|| \delta \bm{\epsilon} ||$',)
plt.ylabel(r'$|| \delta \bm{\sigma}  - \bm{K} : \delta \bm{\epsilon} || / || \delta \bm{\epsilon} ||$')

plt.savefig('consistency_check.svg')
plt.show()
