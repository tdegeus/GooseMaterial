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
#include <cppmat/tensor3.h>

using T2  = cppmat::tensor3_2 <double>;
using T2s = cppmat::tensor3_2s<double>;
using T2d = cppmat::tensor3_2d<double>;
using T4  = cppmat::tensor3_4 <double>;

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

};

// ========================================= IMPLEMENTATION ========================================

LinearElastic::LinearElastic( double K, double G ) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

T2s  LinearElastic::stress(const T2s &eps)
{
  double eps_m,sig_m;
  T2s eps_d,sig_d;
  T2d I;

  // second order identity tensor
  I      = cppmat::identity3_2();

  // decompose strain: hydrostatic part, deviatoric part
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;

  // constitutive response
  sig_m  = 3. * m_K * eps_m;
  sig_d  = 2. * m_G * eps_d;

  // combine volumetric and deviatoric stress
  return sig_m * I + sig_d ;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> LinearElastic::tangent_stress(const T2s &eps)
{
  double eps_m,sig_m;
  T2s eps_d,sig_d,sig;
  T2d I;

  // stress
  // ------

  // second order identity tensor
  I      = cppmat::identity3_2();

  // decompose strain: hydrostatic part, deviatoric part
  eps_m  = eps.trace() / 3.;
  eps_d  = eps - eps_m * I;

  // constitutive response
  sig_m  = 3. * m_K * eps_m;
  sig_d  = 2. * m_G * eps_d;

  // combine volumetric and deviatoric stress
  sig    = sig_m * I + sig_d ;

  // tangent
  // -------

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cppmat::identity3_4d();
  T4 II  = cppmat::identity3_II();

  // initialize tangent as the elasticity tensor
  T4 K4  = m_K * II + 2. * m_G * I4d;

  return std::make_tuple(K4,sig);
}

// -------------------------------------------------------------------------------------------------

}
