/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid

Overview
--------

class ElasticPlasticPotential
|- stress
|- find
|- eps_y

Description
-----------

Elastic plastic material for amorphous solids.

Suggested references
--------------------

*   The code + comments below.
*   docs/ElasticPlasticPotential/ElasticPlasticPotential.pdf

================================================================================================= */

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <tuple>
#include <math.h>
#include <cppmat/tensor.h>

using T2  = cppmat::tensor2 <double>;
using T2s = cppmat::tensor2s<double>;
using T2d = cppmat::tensor2d<double>;
using T4  = cppmat::tensor4 <double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

class ElasticPlasticPotential
{
private:

  double              m_K;      // bulk  modulus
  double              m_G;      // shear modulus
  std::vector<double> m_epsy;      // yield strains
  bool                m_smooth; // set potential smooth or not

public:

  // constructor / destructor
 ~ElasticPlasticPotential(){};
  ElasticPlasticPotential(){};
  ElasticPlasticPotential(double K, double G, const std::vector<double> &epsy, bool smooth=true);

  // compute stress at "eps"
  T2s stress(const T2s &eps);

  // post-process functions
  // - find the index of the current yield strain (below "eps_eq")
  size_t find(double     eps_eq);
  size_t find(const T2s &eps   );
  // - return yield strain "i"
  double eps_y(size_t i);

};

// ========================================= IMPLEMENTATION ========================================

ElasticPlasticPotential::ElasticPlasticPotential(
  double K, double G, const std::vector<double> &epsy, bool smooth )
{
  m_K      = K;
  m_G      = G;
  m_smooth = smooth;

  m_epsy.resize(epsy.size());

  for ( size_t i = 0 ; i < epsy.size() ; ++i )
    m_epsy[i] = epsy[i];
}

// -------------------------------------------------------------------------------------------------

size_t ElasticPlasticPotential::find(double eps_eq)
{
  for ( size_t i = 0 ; i < m_epsy.size()-1 ; ++i )
    if ( eps_eq >= m_epsy[i] && eps_eq < m_epsy[i+1] )
      return i;

  throw std::runtime_error("Insufficient 'eps_y'");
}

// -------------------------------------------------------------------------------------------------

size_t ElasticPlasticPotential::find(const T2s &eps)
{
  double eps_m,eps_eq;
  T2s eps_d;
  T2d I;

  // second order identity tensor
  I      = cppmat::identity2(3);

  // decompose strain: hydrostatic part, deviatoric part
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;

  // equivalent strain
  eps_eq = std::pow( 2./3. * eps_d.ddot(eps_d) , 0.5 );

  // return index of the yield strain
  return find(eps_eq);
}

// -------------------------------------------------------------------------------------------------

double ElasticPlasticPotential::eps_y(size_t i)
{
  return m_epsy[i];
}

// -------------------------------------------------------------------------------------------------

T2s  ElasticPlasticPotential::stress(const T2s &eps)
{
  double eps_m,sig_m,eps_eq;
  T2s eps_d,sig_d,sig;
  T2d I;

  // second order identity tensor
  I      = cppmat::identity2(3);

  // decompose strain: hydrostatic part, deviatoric part
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;

  // equivalent strain
  eps_eq = std::pow( 2./3. * eps_d.ddot(eps_d) , 0.5 );

  // constitutive response - hydrostatic part
  sig_m  = 3. * m_K * eps_m;

  // equivalent strain zero -> zero deviatoric stress
  if ( eps_eq <= 0. ) return sig_m * I;

  // read current yield strains
  size_t i       = find(eps_eq);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // constitutive response - deviatoric part
  if ( ! m_smooth )
    sig_d = 2.*m_G * ( 1. - eps_min/eps_eq ) * eps_d;
  else
    sig_d = 2.*m_G * ( deps_y/M_PI ) * sin ( M_PI/deps_y * ( eps_eq - eps_min ) ) * eps_d / eps_eq;

  // return full strain tensor
  return sig_m * I + sig_d;
}

// -------------------------------------------------------------------------------------------------

}
