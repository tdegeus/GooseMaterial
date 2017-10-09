/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Elasto-visco-plastic material model, whereby the elasticity is based on linear elasticity (i.e. a
linear relationship between the linear strain and the Cauchy stress).

Suggested references
--------------------

*   The code + comments below.
*   docs/Metal/LinearStrain/ElastoViscoPlastic/readme.pdf
*   Former internal code: GooseFEM / mat3202

================================================================================================= */

#include <assert.h>
#include <tuple>
#include <cppmat/tensor3.h>

#warning "GooseMaterial/Metal/LinearStrain/ElastoViscoPlastic/Cartesian3d.h : first usage, careful check then remove this message"

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace ElastoViscoPlastic {
namespace Cartesian3d {

namespace cm = cppmat::cartesian3d;

using T2s = cm::tensor2s<double>;
using T2d = cm::tensor2d<double>;
using T4  = cm::tensor4 <double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:

  double m_K;      // material parameter : bulk  modulus
  double m_G;      // material parameter : shear modulus
  double m_sig0;   // material parameter : 'yield' stress
  double m_gamma0; // material parameter : reference plastic strain rate
  double m_m;      // material parameter : rate exponent
  T2s    m_eps;    // history  parameter : strain tensor
  T2s    m_eps_n;  // history  parameter : strain tensor at last increment
  T2s    m_epse;   // history  parameter : elastic strain tensor
  T2s    m_epse_n; // history  parameter : elastic strain tensor at last increment
  double m_ep;     // history  parameter : accumulated plastic strain
  double m_ep_n;   // history  parameter : accumulated plastic strain at last increment

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, const double dt, bool tangent);

  // compute plastic multiplier for a given trial state and the time-step
  double plastic_multiplier(double sig_eq, double dt);

public:

  // constructor / destructor
 ~Material(){};
  Material(){};
  Material(double K, double G, double sig0, double gamma0, double m=1.);

  // compute stress(+tangent) at "eps", depending on the history stored in this class
  T2s                stress        (const T2s &eps, const double dt);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps, const double dt);

  // update history
  void increment();

};

// ========================================= IMPLEMENTATION ========================================

Material::Material(
  double K, double G, double sig0, double gamma0, double m ) :
  m_K(K), m_G(G), m_sig0(sig0), m_gamma0(gamma0), m_m(m)
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

void Material::increment()
{
  m_eps_n  = m_eps ;
  m_epse_n = m_epse;
  m_ep_n   = m_ep  ;
}

// -------------------------------------------------------------------------------------------------

T2s  Material::stress(const T2s &eps, const double dt)
{
  double epse_m,sig_m,sig_eq,dgamma;
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

  // return map
  if ( sig_eq > 0.0 )
  {
    // - plastic multiplier
    dgamma = plastic_multiplier( sig_eq , dt );
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

std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps, const double dt)
{
  double epse_m,sig_m,sig_eq,dgamma;
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

  // return map
  if ( sig_eq > 0.0 )
  {
    // - plastic multiplier
    dgamma = plastic_multiplier( sig_eq , dt );
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
  if ( sig_eq > 0.0 )
  {
    // - update the tangent (1/2)
    K4 -= 6. * std::pow(m_G,2.) * dgamma/sig_eq * I4d;
    // - update the tangent (2/2)
    K4 += 4. * std::pow(m_G,2.) * (
      dgamma/sig_eq - 1./(3.*m_G + m_m*m_sig0/(m_gamma0*dt) * std::pow(dgamma/(m_gamma0*dt),m_m-1.))
    ) * N.dyadic(N);
  }

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

double Material::plastic_multiplier(double sig_eq, double dt)
{
  // linear stress sensitivity ( m == 1 )
  // ------------------------------------

  if ( std::abs(m_m-1.) < 1.e-6 )
    return sig_eq / ( 3.*m_G + m_sig0/(m_gamma0*dt) );

  // non-linear stress sensitivity ( m != 1 )
  // ----------------------------------------

  int    i      = 0;
  double dgamma = 0.0;
  double R,dR;

  // loop until residual vanishes
  while ( true )
  {

    // - residual
    R  = dgamma - m_gamma0 * dt * std::pow( (sig_eq-3.*m_G*dgamma)/m_sig0 , 1./m_m);

    // - derivative of the residual
    dR = 1. + 3.*m_G*m_gamma0*dt/m_m * std::pow( (sig_eq-3.*m_G*dgamma)/m_sig0 , 1./m_m-1.);

    // - update plastic multiplier
    dgamma -= R/dR;

    // - check for convergence
    if ( std::abs(R/m_gamma0) < 1.0e-6 )
    {
      assert( dgamma >= 0.0 );
      return dgamma;
    }

    // - limit maximum number of iterations
    if ( i>20 )
      throw std::runtime_error("Return-map not succeeded");

    // - update the iteration counter
    ++i;

  }
}

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
