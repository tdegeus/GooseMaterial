/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Elastic plastic material for amorphous solids.

N.B. The elastic part of this model is also implemented in:
"AmorphousSolid/LinearStrain/Elastic/Cartesian3d.h"

Suggested references
--------------------

*   The code + comments below.
*   docs/AmorphousSolid/LinearStrain/ElastoPlastic/readme.pdf

================================================================================================= */

#ifndef GOOSEMATERIAL_AMORPHOUSSOLID_LINEARSTRAIN_ELASTOPLASTIC_SMOOTH_CARTESIAN3D_H
#define GOOSEMATERIAL_AMORPHOUSSOLID_LINEARSTRAIN_ELASTOPLASTIC_SMOOTH_CARTESIAN3D_H

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <tuple>
#include <math.h>
#include <cppmat/tensor3.h>

// -------------------------------------------------------------------------------------------------

namespace GooseMaterial {
namespace AmorphousSolid {
namespace LinearStrain {
namespace ElastoPlastic {
namespace Smooth {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

namespace cm   = cppmat::cartesian3d;
using     T2s  = cm::tensor2s<double>;
using     T2d  = cm::tensor2d<double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:
  double              m_K;    // bulk  modulus
  double              m_G;    // shear modulus
  std::vector<double> m_epsy; // yield strains

public:
  Material(){};
  Material(double K, double G, const std::vector<double> &epsy={}, bool init_elastic=true);

  // compute stress or the energy at "Eps"
  T2s    stress(const T2s &Eps);
  double energy(const T2s &Eps);
  // find the index of the current yield strain, return the yield strain at this index
  size_t find  (double epsd);
  double eps_y (size_t i   );
};

// ===================================== IMPLEMENTATION : CORE =====================================

inline Material::Material(double K, double G, const std::vector<double> &epsy, bool init_elastic)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;

  // copy input - yield strains
  std::vector<double> vec = epsy;
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
  for ( auto &j : vec ) { m_epsy[i] = j; ++i; }

  // check the number of yield strains
  if ( m_epsy.size() < 2 )
    throw std::runtime_error("Specify at least two yield strains 'eps_y'");
}

// -------------------------------------------------------------------------------------------------

inline double Material::eps_y(size_t i)
{
  return m_epsy[i];
}

// -------------------------------------------------------------------------------------------------

inline size_t Material::find(double epsd)
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
    else                     r = i;

    // set new search index
    i = ( r + l ) / 2;
    i = std::max(i,z);
  }
}

// -------------------------------------------------------------------------------------------------

inline T2s Material::stress(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/3.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // no deviatoric strain -> return stress tensor with only hydrostatic part
  if ( epsd <= 0. ) return (m_K*epsm) * I;

  // read current yield strains
  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // return full strain tensor
  return (m_K*epsm)*I + ((m_G/epsd)*(deps_y/M_PI)*sin(M_PI/deps_y*(epsd-eps_min)))*Epsd;
}

// -------------------------------------------------------------------------------------------------

inline double Material::energy(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/3.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // energy for the hydrostatic part
  double U = 3./2. * m_K * std::pow(epsm,2.);

  // read current yield strain
  size_t i       = find(epsd);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // deviatoric part of the energy
  double V = -2.*m_G * std::pow(deps_y/M_PI,2.) * ( 1. + cos( M_PI/deps_y * (epsd-eps_min) ) );

  // return total energy
  return U + V;
}

// =================================================================================================

}}}}}} // namespace ...

#endif
