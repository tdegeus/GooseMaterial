/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMesh

Overview
--------

class PlasticLinearElastic
|- stress
|- stress_tangent
|- tangent
|- increment

Description
-----------

Elasto-plastic material model, whereby the elasticity is based on linear elasticity (i.e. a linear
relationship between the linear strain and the Cauchy stress).

Suggested references
--------------------

*   The code + comments below.
*   docs/PlasticLinearElastic/PlasticLinearElastic.pdf
*   Former internal code: GooseFEM / mat2002

================================================================================================= */

#include <tuple>
#include <cppmat/tensor.h>

using T2  = cppmat::tensor2 <double>;
using T2s = cppmat::tensor2s<double>;
using T2d = cppmat::tensor2d<double>;
using T4  = cppmat::tensor4 <double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

class PlasticLinearElastic
{
private:
  double m_K      ; // material parameter : bulk  modulus
  double m_G      ; // material parameter : shear modulus
  double m_sigy0  ; // material parameter : initial yield stress
  double m_H      ; // material parameter : hardening modulus
  double m_m      ; // material parameter : hardening exponent
  T2s    m_eps    ; // history  parameter : strain tensor
  T2s    m_eps_n  ; // history  parameter : strain tensor at last increment
  T2s    m_epse   ; // history  parameter : elastic strain tensor
  T2s    m_epse_n ; // history  parameter : elastic strain tensor at last increment
  double m_ep     ; // history  parameter : accumulated plastic strain
  double m_ep_n   ; // history  parameter : accumulated plastic strain at last increment

public:
  // constructor / destructor
 ~PlasticLinearElastic(){};
  PlasticLinearElastic(){};
  PlasticLinearElastic(double K, double G, double sigy0, double H, double m=1.);

  // compute stress(+tangent) at "eps", depending on the history stored in this class
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> stress_tangent(const T2s &eps);
  T4                 tangent       (const T2s &eps);

  // update history
  void increment();

  // perform actual computations
  std::tuple<T4,T2s> f_compute(const T2s &eps, bool stress_only=false);
};

// ========================================= IMPLEMENTATION ========================================

PlasticLinearElastic::PlasticLinearElastic(
  double K, double G, double sigy0, double H, double m ) :
  m_K(K), m_G(G), m_sigy0(sigy0), m_H(H), m_m(m)
{
  // resize history tensors
  m_eps   .resize(3);
  m_eps_n .resize(3);
  m_epse  .resize(3);
  m_epse_n.resize(3);

  // initialize stress/strain free state
  m_eps_n .zeros();
  m_epse_n.zeros();
  m_ep_n = 0.0;
}

// -------------------------------------------------------------------------------------------------

void PlasticLinearElastic::increment()
{
  m_eps_n  = m_eps ;
  m_epse_n = m_epse;
  m_ep_n   = m_ep  ;
}

// -------------------------------------------------------------------------------------------------

T2s  PlasticLinearElastic::stress(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = this->f_compute(eps,true);

  return sig;
}

// -------------------------------------------------------------------------------------------------

T4   PlasticLinearElastic::tangent(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = this->f_compute(eps,false);

  return K4;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> PlasticLinearElastic::stress_tangent(const T2s &eps)
{
  return this->f_compute(eps,false);
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> PlasticLinearElastic::f_compute(const T2s &eps, bool stress_only)
{
  double dgamma = 0.0;
  double epse_m,sig_m,sig_eq,dH,phi;
  T2s epse_d,sig_d,sig,N;
  T2d I;

  // stress response
  // ---------------

  // second order identity tensor
  I      = cppmat::identity2(3);

  // total strain, trial elastic strain, and accumulated plastic strain
  m_eps  = eps;
  m_epse = m_epse_n + ( eps - m_eps_n );
  m_ep   = m_ep_n;

  // decompose trial elastic strain: hydrostatic part, deviatoric part
  epse_m = m_epse.trace()/3.;
  epse_d = m_epse - epse_m * I;

  // trial stress: hydrostatic and deviatoric part, and equivalent stress
  sig_m  = 3. * m_K * epse_m;
  sig_d  = 2. * m_G * epse_d;
  sig_eq = std::pow( 1.5 * sig_d.ddot(sig_d) , 0.5 );

  // yield surface, initialize return-map variables
  phi    = sig_eq - m_sigy0 - m_H * std::pow( m_ep_n , m_m );
  dH     = m_H;
  dgamma = 0.0;

  // return map
  if ( phi > 0.0 )
  {
    // - determine plastic multiplier
    //   (a) linear hardening
    if ( std::abs(m_m-1.) < 1.e-6 )
    {
      dgamma = phi / ( 3.*m_G + m_H );
    }
    //   (b) non-linear hardening
    else
    {
      int    i = 0;
      double d,R;
      // -- loop until residual vanishes: do { ... } while ( residual > norm )
      do
      {
        // --- residual
        R = sig_eq - 3.*m_G*dgamma - m_sigy0 - m_H * std::pow( m_ep_n+dgamma , m_m );
        // --- hardening slope (avoid zero-division)
        if ( std::abs( m_ep_n+dgamma ) < 1.e-6 ) {
          d  = R / ( -3.*m_G );
        }
        else {
          dH = m_m * m_H * std::pow( m_ep_n+dgamma , m_m-1. );
          d  = R / ( -3.*m_G - dH );
        }
        // --- avoid zero devision
        if ( m_ep_n+dgamma-d <= 0. )
          d = ( m_ep_n+dgamma )*.9;
        // --- update plastic multiplier
        dgamma -= d;
        // --- limit maximum number of iterations
        if ( i>20 ) throw std::runtime_error("Return-map not succeeded");
        ++i;
      }
      while ( std::abs(R/m_sigy0) > 1.e-6 );
      // -- double check physical admissibility
      if ( dgamma < 0.0 ) throw std::runtime_error("Negative plastic multiplier found");
    }

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

  if ( stress_only ) {
    T4 K4(0);
    return std::make_tuple(K4,sig);
  }

  // fourth order deviator unit tensor (i.e. A_d = I4d : A)
  T4 I4d = cppmat::identity4d (3);
  T4 II  = cppmat::identity4II(3);

  // initialize tangent as the elasticity tensor
  T4 K4  = m_K * II + 2. * m_G * I4d;

  // plastic part: only when yielding
  if ( phi > 0.0 )
  {
    // - tangential hardening modulus
    dH  = m_m * m_H * std::pow( m_ep , m_m-1. );
    // - update the tangent (1/2)
    K4 -= 6. * std::pow(m_G,2.) * dgamma/sig_eq * I4d;
    // - update the tangent (2/2)
    K4 += 4. * std::pow(m_G,2.) * ( dgamma/sig_eq - 1./(3.*m_G+dH) ) * N.dyadic(N);
  }

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

}
