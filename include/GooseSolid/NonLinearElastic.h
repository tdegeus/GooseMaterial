/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid

Overview
--------

class NonLinearElastic
|- stress
|- tangent_stress
|- tangent

Description
-----------

Non-linear elasticity: a specific non-linear relationship between the linear strain and the Cauchy
stress.

Suggested references
--------------------

*   The code + comments below.
*   docs/NonLinearElastic/NonLinearElastic.pdf
*   Former internal code: GooseFEM / mat1101

================================================================================================= */

#include <tuple>
#include <cppmat/tensor.h>

using T2  = cppmat::tensor2 <double>;
using T2s = cppmat::tensor2s<double>;
using T2d = cppmat::tensor2d<double>;
using T4  = cppmat::tensor4 <double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

class NonLinearElastic
{
private:

  double m_K;    // material parameter : bulk modulus
  double m_sig0; // material parameter : reference stress
  double m_eps0; // material parameter : reference strain
  double m_n;    // material parameter : exponent

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

public:

  // constructor / destructor
 ~NonLinearElastic(){};
  NonLinearElastic(){};
  NonLinearElastic(double K, double sig0, double eps0, double m=1.);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);
  T4                 tangent       (const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

NonLinearElastic::NonLinearElastic( double K, double sig0, double eps0, double m ) :
  m_K(K), m_sig0(sig0), m_eps0(eps0), m_n(m)
{
}

// -------------------------------------------------------------------------------------------------

T2s  NonLinearElastic::stress(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = compute(eps,false);

  return sig;
}

// -------------------------------------------------------------------------------------------------

T4   NonLinearElastic::tangent(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = compute(eps,true);

  return K4;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> NonLinearElastic::tangent_stress(const T2s &eps)
{
  return compute(eps,true);
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> NonLinearElastic::compute(const T2s &eps, bool tangent)
{
  double eps_m,sig_m,eps_eq;
  T2s eps_d,sig;
  T2d I;
  T2s sig_d(3, 0.0);

  // stress
  // ------

  // second order identity tensor
  I      = cppmat::identity2(3);

  // decompose strain: hydrostatic part, deviatoric part, equivalent strain
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;
  eps_eq = std::pow( 2./3. * eps_d.ddot(eps_d) , 0.5 );

  // hydrostatic stress
  sig_m  = 3. * m_K * eps_m;

  // deviatoric stress
  if ( eps_eq != 0.0 )
    sig_d = 2./3. * m_sig0/std::pow(m_eps0,m_n) * std::pow(eps_eq,m_n-1.) * eps_d;

  // combine volumetric and deviatoric stress
  sig = sig_m * I + sig_d ;

  // tangent
  // -------

  // compute only stress: allocate empty tangent (without any element), and return
  if ( ! tangent ) {
    T4 K4(0);
    return std::make_tuple(K4,sig);
  }

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cppmat::identity4d (3);
  T4 II  = cppmat::identity4II(3);

  // hydrostatic part
  T4 K4  = m_K * II;

  // deviator part
  if ( eps_eq != 0.0 )
    K4 += 2./3. * m_sig0/std::pow(m_eps0,m_n)
       * (2./3.*(m_n-1.)*std::pow(eps_eq,m_n-3.)*eps_d.dyadic(eps_d) + std::pow(eps_eq,m_n-1.)*I4d);
  else
    K4 += I4d;

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

}
