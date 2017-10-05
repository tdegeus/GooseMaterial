/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Elastic plastic material for amorphous solids.

Suggested references
--------------------

*   The code + comments below.
*   docs/AmorphousSolid/LinearStrain/ElastoPlastic/readme.pdf

================================================================================================= */

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <tuple>
#include <math.h>
#include <cppmat/tensor3.h>

namespace GooseMaterial {
namespace AmorphousSolid {
namespace LinearStrain {
namespace ElastoPlastic {
namespace Cartesian3d {

using T2s = cppmat::tensor3_2s<double>;
using T2d = cppmat::tensor3_2d<double>;

T2d    I    = cppmat::identity3_2();
double ndim = 3.;

// ============================================ OVERVIEW ===========================================

class Material
{
private:

  double              m_K;        // bulk  modulus
  double              m_G;        // shear modulus
  std::vector<double> m_epsy;     // yield strains
  bool                m_elastic;  // material is elastic or not

public:

  // constructor / destructor
 ~Material(){};
  Material(){};
  Material(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // compute stress at "Eps"
  T2s stress(const T2s &Eps);

  // post-process functions
  // - return is material is elastic or not
  bool   elastic();
  // - find the index of the current yield strain (below "eps_d")
  size_t find(double     epsd);
  size_t find(const T2s &Eps );
  // - return yield strain "i"
  double eps_y(size_t i);
  // - return (hydrostatic/deviatoric) equivalent stress/strain
  double eps_eq(const T2s &Eps);
  double eps_m (const T2s &Eps);
  double eps_d (const T2s &Eps);
  double sig_eq(const T2s &Sig);
  double sig_m (const T2s &Sig);
  double sig_d (const T2s &Sig);
  // - return the strain energy (or its hydrostatic or deviatoric component)
  double energy  (double epsm, double epsd);
  double energy_m(double epsm);
  double energy_d(double epsd);
  double energy  (const T2s &Eps);
  double energy_m(const T2s &Eps);
  double energy_d(const T2s &Eps);

};

// ===================================== IMPLEMENTATION : CORE =====================================

Material::Material(double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;

  // copy input - yield strains
  std::vector<double> vec = epsy;
  // sort input
  std::sort( vec.begin() , vec.end() );

  // elastic : no yield strains
  if ( vec.size() == 0 ) {
    m_elastic = true;
    return;
  }

  // set plastic
  m_elastic = false;

  // check to add item to have an initial elastic response
  if ( init_elastic )
    if ( vec[0] == -vec[1] )
      init_elastic = false;

  // copy input
  // - counters
  size_t N = vec.size();
  size_t i = 0;
  // - add yield strain to have an initial elastic response
  if ( init_elastic ) ++N;
  // - allocate
  m_epsy.resize(N);
  // - add yield strain to have an initial elastic response
  if ( init_elastic ) { m_epsy[i] = -vec[0]; ++i; }
  // - copy the rest
  for ( auto &j : vec ) {
    m_epsy[i] = j; ++i;
  }

  // check the number of yield strains
  if ( m_epsy.size() < 2 )
    throw std::runtime_error("Specify at least two yield strains 'eps_y'");
}

// -------------------------------------------------------------------------------------------------

size_t Material::find(double epsd)
{
  // check extremes
  if ( epsd < m_epsy.front() or epsd >= m_epsy.back() )
    throw std::runtime_error("Insufficient 'eps_y'");

  // set initial search bounds and index
  size_t n = m_epsy.size()-1;  // upper-bound
  size_t z = 1;                // lower-bound
  size_t l = 0;                // left-bound
  size_t r = n;                // right-bound
  size_t i = r / 2;            // average

  // loop until found
  while ( true )
  {
    // check if found, unroll once to speed-up
    if ( epsd >= m_epsy[i-1] and epsd < m_epsy[i  ] ) return i-1;
    if ( epsd >= m_epsy[i  ] and epsd < m_epsy[i+1] ) return i;
    if ( epsd >= m_epsy[i+1] and epsd < m_epsy[i+2] ) return i+1;

    // correct the left- and right-bound
    if ( epsd >= m_epsy[i] ) l = i;
    else                      r = i;

    // set new search index
    i = ( r + l ) / 2;
    i = std::max(i,z);
  }
}

// -------------------------------------------------------------------------------------------------

T2s Material::stress(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  double epsm  = Eps.trace()/ndim;
  T2s    Epsd  = Eps - epsm*I;
  double epsd = std::pow( .5*Epsd.ddot(Epsd) , 0.5 );

  // elastic: return full stress tensor
  if ( m_elastic  ) return (m_K*epsm) * I + m_G * Epsd;

  // zero equivalent strain -> return stress tensor with only hydrostatic part
  if ( epsd <= 0. ) return (m_K*epsm) * I;

  // read current yield strains
  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // return full strain tensor
  return (m_K*epsm)*I + ((m_G/epsd)*(deps_y/M_PI)*sin(M_PI/deps_y*(epsd-eps_min)))*Epsd;
}

// ================================= IMPLEMENTATION : POST-PROCESS =================================

bool Material::elastic()
{
  return m_elastic;
}

// -------------------------------------------------------------------------------------------------

size_t Material::find(const T2s &Eps)
{
  return find(eps_d(Eps));
}

// -------------------------------------------------------------------------------------------------

double Material::eps_y(size_t i)
{
  return m_epsy[i];
}

// -------------------------------------------------------------------------------------------------

double Material::eps_eq(const T2s &Eps)
{
  return std::pow( .5*Eps.ddot(Eps) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::eps_m(const T2s &Eps)
{
  return Eps.trace()/ndim;
}

// -------------------------------------------------------------------------------------------------

double Material::eps_d(const T2s &Eps)
{
  double epsm = Eps.trace()/ndim;
  T2s    Epsd = Eps - epsm*I;

  return std::pow( .5*Epsd.ddot(Epsd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_eq(const T2s &Sig)
{
  return std::pow( .5*Sig.ddot(Sig) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_m(const T2s &Sig)
{
  return Sig.trace()/ndim;
}

// -------------------------------------------------------------------------------------------------

double Material::sig_d(const T2s &Sig)
{
  double sigm = Sig.trace()/ndim;
  T2s    Sigd = Sig - sigm*I;

  return std::pow( .5*Sigd.ddot(Sigd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(double epsm)
{
  return ndim/2. * m_K * std::pow( epsm , 2. );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(const T2s &Eps)
{
  return energy_m(eps_m(Eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy_d(double epsd)
{
  if ( m_elastic  ) return m_G * std::pow( epsd , 2. );
  if ( epsd <= 0. ) return 0.0;

  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  return -2. * m_G * std::pow( deps_y/M_PI , 2. ) * ( 1. + cos( M_PI/deps_y * (epsd-eps_min) ) );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_d(const T2s &Eps)
{
  return energy_d(eps_d(Eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy(double epsm, double epsd)
{
  return energy_m(epsm) + energy_d(epsd);
}

// -------------------------------------------------------------------------------------------------

double Material::energy(const T2s &Eps)
{
  return energy_m(eps_m(Eps)) + energy_d(eps_d(Eps));
}

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
