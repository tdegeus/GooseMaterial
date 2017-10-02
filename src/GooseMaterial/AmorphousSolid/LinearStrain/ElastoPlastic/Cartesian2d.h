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
#include <cppmat/tensor2.h>

#warning "GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h : first usage, careful check then remove this message"

namespace GooseMaterial {
namespace AmorphousSolid {
namespace LinearStrain {
namespace ElastoPlastic {
namespace Cartesian2d {

using T2s = cppmat::tensor2_2s<double>;
using T2d = cppmat::tensor2_2d<double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:

  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:

  // constructor / destructor
 ~Material(){};
  Material(){};
  Material(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // compute stress at "eps"
  T2s stress(const T2s &eps);

  // post-process functions
  // - find the index of the current yield strain (below "eps_eq")
  size_t find(double     epseq);
  size_t find(const T2s &eps  );
  // - return yield strain "i"
  double eps_y(size_t i);
  // - return hydrostatic/equivalent stress/strain
  double eps_eq(const T2s &eps);
  double eps_m (const T2s &eps);
  double sig_eq(const T2s &sig);
  double sig_m (const T2s &sig);
  // - return the strain energy (or its hydrostatic or deviatoric component)
  double energy   (double     epsm, double epseq);
  double energy_m (double     epsm              );
  double energy_eq(double     epseq             );
  double energy   (const T2s &eps               );
  double energy_m (const T2s &eps               );
  double energy_eq(const T2s &eps               );

};

// ===================================== IMPLEMENTATION : CORE =====================================

Material::Material(double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;

  // copy input - yield strains
  std::vector<double> vec = epsy;
  // initialize: no yield strain -> only elastic
  if ( vec.size() == 0 ) vec = { 1000000. };
  // sort input
  std::sort( vec.begin() , vec.end() );

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

size_t Material::find(double epseq)
{
  // check extremes
  if ( epseq < m_epsy.front() or epseq >= m_epsy.back() )
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
    if ( epseq >= m_epsy[i-1] and epseq < m_epsy[i  ] ) return i-1;
    if ( epseq >= m_epsy[i  ] and epseq < m_epsy[i+1] ) return i;
    if ( epseq >= m_epsy[i+1] and epseq < m_epsy[i+2] ) return i+1;

    // correct the left- and right-bound
    if ( epseq >= m_epsy[i] ) l = i;
    else                      r = i;

    // set new search index
    i = ( r + l ) / 2;
    i = std::max(i,z);
  }
}

// -------------------------------------------------------------------------------------------------

T2s Material::stress(const T2s &eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I     = cppmat::identity2_2();
  double epsm  = eps.trace()/2.;
  T2s    epsd  = eps - epsm*I;
  double epseq = std::pow( .5*epsd.ddot(epsd) , 0.5 );

  // constitutive response - hydrostatic part
  double sigm  = m_K * epsm;

  // equivalent strain zero -> zero deviatoric stress
  if ( epseq <= 0. ) return sigm * I;

  // read current yield strains
  size_t i       = find(epseq);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // constitutive response - deviatoric part
  T2s sigd = ( ( m_G/epseq ) * ( deps_y/M_PI ) * sin ( M_PI/deps_y * (epseq-eps_min) ) ) * epsd;

  // return full strain tensor
  return sigm*I + sigd;
}

// ================================= IMPLEMENTATION : POST-PROCESS =================================

size_t Material::find(const T2s &eps)
{
  return find(eps_eq(eps));
}

// -------------------------------------------------------------------------------------------------

double Material::eps_y(size_t i)
{
  return m_epsy[i];
}

// -------------------------------------------------------------------------------------------------

double Material::eps_eq(const T2s &eps)
{
  T2d    I    = cppmat::identity2_2();
  double epsm = eps.trace()/2.;
  T2s    epsd = eps - epsm*I;

  return std::pow( .5*epsd.ddot(epsd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::eps_m(const T2s &eps)
{
  return eps.trace()/2.;
}

// -------------------------------------------------------------------------------------------------

double Material::sig_eq(const T2s &sig)
{
  T2d    I    = cppmat::identity2_2();
  double sigm = sig.trace()/2.;
  T2s    sigd = sig - sigm*I;

  return std::pow( .5*sigd.ddot(sigd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_m(const T2s &sig)
{
  return sig.trace()/2.;
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(double epsm)
{
  return m_K * std::pow( epsm , 2. );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(const T2s &eps)
{
  return energy_m(eps_m(eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy_eq(double epseq)
{
  size_t i       = find(epseq);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  return -2. * m_G * std::pow( deps_y/M_PI , 2. ) * ( 1. + cos( M_PI/deps_y * (epseq-eps_min) ) );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_eq(const T2s &eps)
{
  return energy_eq(eps_eq(eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy(double epsm, double epseq)
{
  return energy_m(epsm) + energy_eq(epseq);
}

// -------------------------------------------------------------------------------------------------

double Material::energy(const T2s &eps)
{
  return energy_m(eps_m(eps)) + energy_eq(eps_eq(eps));
}

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
