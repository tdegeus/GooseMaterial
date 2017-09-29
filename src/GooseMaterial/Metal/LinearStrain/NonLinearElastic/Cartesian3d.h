/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Non-linear elasticity: a specific non-linear relationship between the linear strain and the Cauchy
stress.

Suggested references
--------------------

*   The code + comments below.
*   docs/Metal/LinearStrain/NonLinearElastic/readme.pdf
*   Former internal code: GooseFEM / mat1101

================================================================================================= */

#include <tuple>
#include <cppmat/tensor3.h>

#warning "GooseMaterial/Metal/LinearStrain/NonLinearElastic/Cartesian3d.h : first usage, careful check then remove this message"

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace NonLinearElastic {

using T2  = cppmat::tensor3_2 <double>;
using T2s = cppmat::tensor3_2s<double>;
using T2d = cppmat::tensor3_2d<double>;
using T4  = cppmat::tensor3_4 <double>;

// ============================================ OVERVIEW ===========================================

class Material
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
 ~Material(){};
  Material(){};
  Material(double K, double sig0, double eps0, double m=1.);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

Material::Material( double K, double sig0, double eps0, double m ) :
  m_K(K), m_sig0(sig0), m_eps0(eps0), m_n(m)
{
}

// -------------------------------------------------------------------------------------------------

T2s  Material::stress(const T2s &eps)
{
  double eps_m,sig_m,eps_eq;
  T2s eps_d;
  T2d I;
  T2s sig_d(0.0);

  // second order identity tensor
  I      = cppmat::identity3_2();

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
  return sig_m * I + sig_d ;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps)
{
  double eps_m,sig_m,eps_eq;
  T2s eps_d,sig;
  T2d I;
  T2s sig_d(0.0);

  // stress
  // ------

  // second order identity tensor
  I      = cppmat::identity3_2();

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

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cppmat::identity3_4d();
  T4 II  = cppmat::identity3_II();

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

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
