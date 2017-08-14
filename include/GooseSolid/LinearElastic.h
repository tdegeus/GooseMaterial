/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid

Overview
--------

class LinearElastic
|- stress
|- tangent_stress
|- tangent

Description
-----------

Linear elasticity: a linear relationship between the linear strain and the Cauchy stress.

Suggested references
--------------------

*   The code + comments below.
*   docs/LinearElastic/LinearElastic.pdf
*   Former internal code: GooseFEM / mat1001

================================================================================================= */

#include <tuple>
#include <cppmat/tensor.h>

using T2  = cppmat::tensor2 <double>;
using T2s = cppmat::tensor2s<double>;
using T2d = cppmat::tensor2d<double>;
using T4  = cppmat::tensor4 <double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

class LinearElastic
{
private:

  double m_K; // material parameter : bulk  modulus
  double m_G; // material parameter : shear modulus

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

public:

  // constructor / destructor
 ~LinearElastic(){};
  LinearElastic(){};
  LinearElastic(double K, double G);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);
  T4                 tangent       (const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

LinearElastic::LinearElastic( double K, double G ) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

T2s  LinearElastic::stress(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = compute(eps,false);

  return sig;
}

// -------------------------------------------------------------------------------------------------

T4   LinearElastic::tangent(const T2s &eps)
{
  T2s sig;
  T4  K4;

  std::tie(K4,sig) = compute(eps,true);

  return K4;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> LinearElastic::tangent_stress(const T2s &eps)
{
  return compute(eps,true);
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> LinearElastic::compute(const T2s &eps, bool tangent)
{
  double eps_m,sig_m;
  T2s eps_d,sig_d,sig;
  T2d I;

  // stress
  // ------

  // second order identity tensor
  I      = cppmat::identity2(3);

  // decompose strain: hydrostatic part, deviatoric part
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;

  // constitutive response
  sig_m  = 3. * m_K * eps_m;
  sig_d  = 2. * m_G * eps_d;

  // combine volumetric and deviatoric stress
  sig    = sig_m  * I + sig_d ;

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

  // initialize tangent as the elasticity tensor
  T4 K4  = m_K * II + 2. * m_G * I4d;

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

}
