/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Elasto-plastic material model with power-law hardening of the initial yield stress. The elasticity
is based on linear elasticity (i.e. a linear relationship between the linear strain and the Cauchy
stress).

Suggested references
--------------------

*   The code + comments below.
*   docs/Metal/LinearStrain/ElastoPlastic/readme.pdf
*   Former internal code: GooseFEM / mat2002

================================================================================================= */

#ifndef GOOSEMATERIAL_METAL_LINEARSTRAIN_ELASTOPLASTIC_CARTESIAN3D_H
#define GOOSEMATERIAL_METAL_LINEARSTRAIN_ELASTOPLASTIC_CARTESIAN3D_H

#include <assert.h>
#include <tuple>
#include <cppmat/tensor3.h>

// -------------------------------------------------------------------------------------------------

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace ElastoPlastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

namespace cm  = cppmat::cartesian3d;
using     T2s = cm::tensor2s<double>;
using     T2d = cm::tensor2d<double>;
using     T4  = cm::tensor4 <double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:
  double m_K;      // material parameter : bulk  modulus
  double m_G;      // material parameter : shear modulus
  double m_sigy0;  // material parameter : initial yield stress
  double m_H;      // material parameter : hardening modulus
  double m_m;      // material parameter : hardening exponent
  T2s    m_eps;    // history  parameter : strain tensor
  T2s    m_eps_n;  // history  parameter : strain tensor at last increment
  T2s    m_epse;   // history  parameter : elastic strain tensor
  T2s    m_epse_n; // history  parameter : elastic strain tensor at last increment
  double m_ep;     // history  parameter : accumulated plastic strain
  double m_ep_n;   // history  parameter : accumulated plastic strain at last increment

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

  // compute the plastic multiplier and tangent hardening modulus for a given trial state
  std::tuple<double,double> plastic_multiplier(double phi, double sig_eq);

public:
  Material(){};
  Material(double K, double G, double sigy0, double H, double m=1.);

  // compute stress(+tangent) at "eps", depending on the history stored in this class
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);

  // update history
  void increment();

};

// ========================================= IMPLEMENTATION ========================================

inline Material::Material(
  double K, double G, double sigy0, double H, double m ) :
  m_K(K), m_G(G), m_sigy0(sigy0), m_H(H), m_m(m)
{
  // initialize stress/strain free state
  m_eps   .zeros();
  m_eps_n .zeros();
  m_epse  .zeros();
  m_epse_n.zeros();
  m_ep   = 0.0;
  m_ep_n = 0.0;
}

// -------------------------------------------------------------------------------------------------

inline void Material::increment()
{
  m_eps_n  = m_eps ;
  m_epse_n = m_epse;
  m_ep_n   = m_ep  ;
}

// -------------------------------------------------------------------------------------------------

inline T2s Material::stress(const T2s &eps)
{
  double epse_m,sig_m,sig_eq,phi,dgamma,dH;
  T2s epse_d,sig_d;
  T2d I;

  // second order identity tensor
  I      = cm::identity2();

  // trial strain: copy total strain, set trial elastic strain and trial accumulated plastic strain
  m_eps  = eps;
  m_epse = m_epse_n + ( eps - m_eps_n );
  m_ep   = m_ep_n;

  // decompose trial elastic strain: hydrostatic part, deviatoric part
  epse_m = m_epse.trace() / 3.;
  epse_d = m_epse - epse_m * I;

  // trial stress: hydrostatic and deviatoric part, equivalent stress
  sig_m  = 3. * m_K * epse_m;
  sig_d  = 2. * m_G * epse_d;
  sig_eq = std::pow( 1.5 * sig_d.ddot(sig_d) , 0.5 );

  // yield surface
  phi    = sig_eq - ( m_sigy0 + m_H * std::pow( m_ep_n , m_m ) );

  // return map
  if ( phi > 0.0 )
  {
    // - plastic multiplier (and tangential hardening modulus)
    std::tie( dgamma , dH ) = plastic_multiplier( phi , sig_eq );
    // - correct trial state
    sig_d  = ( 1.  -  3.*m_G*dgamma/sig_eq ) * sig_d;
    epse_d = sig_d / (2.*m_G);
    m_ep  += dgamma;
  }

  // combine volumetric and deviatoric strain
  m_epse = epse_m * I + epse_d;

  // combine volumetric and deviatoric stress
  return sig_m * I + sig_d ;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps)
{
  double epse_m,sig_m,sig_eq,phi,dgamma,dH;
  T2s epse_d,sig_d,sig,N;
  T2d I;

  // stress
  // ------

  // second order identity tensor
  I      = cm::identity2();

  // trial strain: copy total strain, set trial elastic strain and trial accumulated plastic strain
  m_eps  = eps;
  m_epse = m_epse_n + ( eps - m_eps_n );
  m_ep   = m_ep_n;

  // decompose trial elastic strain: hydrostatic part, deviatoric part
  epse_m = m_epse.trace() / 3.;
  epse_d = m_epse - epse_m * I;

  // trial stress: hydrostatic and deviatoric part, equivalent stress
  sig_m  = 3. * m_K * epse_m;
  sig_d  = 2. * m_G * epse_d;
  sig_eq = std::pow( 1.5 * sig_d.ddot(sig_d) , 0.5 );

  // yield surface
  phi    = sig_eq - ( m_sigy0 + m_H * std::pow( m_ep_n , m_m ) );

  // return map
  if ( phi > 0.0 )
  {
    // - plastic multiplier (and tangential hardening modulus, for tangent)
    std::tie( dgamma , dH ) = plastic_multiplier( phi , sig_eq );
    // - yield surface normal (for tangent)
    N      = 1.5 * sig_d/sig_eq;
    // - correct trial state
    sig_d  = ( 1.  -  3.*m_G*dgamma/sig_eq ) * sig_d;
    epse_d = sig_d / (2.*m_G);
    m_ep  += dgamma;
  }

  // combine volumetric and deviatoric stress/strain
  sig    = sig_m  * I + sig_d ;
  m_epse = epse_m * I + epse_d;

  // tangent
  // -------

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cm::identity4d();
  T4 II  = cm::identityII();

  // initialize tangent as the elasticity tensor
  T4 K4  = m_K * II + 2. * m_G * I4d;

  // plastic part: only when yielding
  if ( phi > 0.0 )
  {
    // - update the tangent (1/2)
    K4 -= 6. * std::pow(m_G,2.) *   dgamma/sig_eq * I4d;
    // - update the tangent (2/2)
    K4 += 4. * std::pow(m_G,2.) * ( dgamma/sig_eq - 1./(3.*m_G+dH) ) * N.dyadic(N);
  }

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<double,double> Material::plastic_multiplier(double phi, double sig_eq)
{
  // linear hardening ( m == 1 )
  // ---------------------------

  if ( std::abs(m_m-1.) < 1.e-6 )
    return std::make_tuple( phi / ( 3.*m_G + m_H ) , m_H );

  // non-linear hardening ( m != 1 )
  // -------------------------------

  int    i      = 0;
  double dgamma = 0.0;
  double d,R,dH;

  // loop until residual vanishes
  while ( true )
  {

    // - residual
    R = sig_eq - 3.*m_G*dgamma - m_sigy0 - m_H * std::pow( m_ep_n+dgamma , m_m );

    // - hardening slope (avoid zero devision, since often "m_m - 1.0 < 0.0")
    if ( std::abs( m_ep_n+dgamma ) < 1.e-8 )
      dH = 0.0;
    else
      dH = m_m * m_H * std::pow( m_ep_n+dgamma , m_m-1. );

    // - actual iterative update of dgamma
    d  = R / ( -3.*m_G - dH );

    // - numerical fix: avoid inadmissible state "ep < 0.0"
    //   (this could happen is the hardening is very non-linear)
    if ( m_ep_n+dgamma-d <= 0. )
      d = ( m_ep_n+dgamma )*.9;

    // - update plastic multiplier
    dgamma -= d;

    // - check for convergence
    if ( std::abs(R/m_sigy0) < 1.0e-6 )
    {
      assert( dgamma >= 0.0 );
      return std::make_tuple( dgamma , dH );
    }

    // - limit maximum number of iterations
    if ( i>20 )
      throw std::runtime_error("Return-map not succeeded");

    // - update the iteration counter
    ++i;

  }
}

// =================================================================================================

}}}}} // namespace ...

#endif
