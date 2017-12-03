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

#ifndef GOOSEMATERIAL_METAL_LINEARSTRAIN_NONLINEARELASTIC_CARTESIAN3D_H
#define GOOSEMATERIAL_METAL_LINEARSTRAIN_NONLINEARELASTIC_CARTESIAN3D_H

#include <tuple>
#include <cppmat/cppmat.h>

// -------------------------------------------------------------------------------------------------

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace NonLinearElastic {
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
  double m_K;    // material parameter : bulk modulus
  double m_sig0; // material parameter : reference stress
  double m_eps0; // material parameter : reference strain
  double m_n;    // material parameter : exponent

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

public:
  Material(){};
  Material(double K, double sig0, double eps0, double m=1.);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

inline Material::Material( double K, double sig0, double eps0, double m ) :
  m_K(K), m_sig0(sig0), m_eps0(eps0), m_n(m)
{
}

// -------------------------------------------------------------------------------------------------

inline T2s  Material::stress(const T2s &eps)
{
  double epsm,sigm,epseq;
  T2s epsd;
  T2d I;
  T2s sigd(0.0);

  // second order identity tensor
  I     = cm::identity2<double>();

  // decompose strain: hydrostatic part, deviatoric part, equivalent strain
  epsm  = eps.trace() / 3.;
  epsd  = eps - epsm * I;
  epseq = std::pow( 2./3. * epsd.ddot(epsd) , 0.5 );

  // hydrostatic stress
  sigm  = ( 3. * m_K ) * epsm;

  // deviatoric stress
  if ( epseq != 0.0 )
    sigd = ( 2./3. * m_sig0/std::pow(m_eps0,m_n) * std::pow(epseq,m_n-1.) ) * epsd;

  // combine volumetric and deviatoric stress
  return sigm * I + sigd ;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps)
{
  double epsm,sigm,epseq;
  T2s epsd,sig;
  T2d I;
  T2s sigd(0.0);

  // stress
  // ------

  // second order identity tensor
  I     = cm::identity2<double>();

  // decompose strain: hydrostatic part, deviatoric part, equivalent strain
  epsm  = eps.trace() / 3.;
  epsd  = eps - epsm * I;
  epseq = std::pow( 2./3. * epsd.ddot(epsd) , 0.5 );

  // hydrostatic stress
  sigm  = ( 3. * m_K ) * epsm;

  // deviatoric stress
  if ( epseq != 0.0 )
    sigd = ( 2./3. * m_sig0/std::pow(m_eps0,m_n) * std::pow(epseq,m_n-1.) ) * epsd;

  // combine volumetric and deviatoric stress
  sig = sigm * I + sigd ;

  // tangent
  // -------

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cm::identity4d<double>();
  T4 II  = cm::identity4II<double>();

  // hydrostatic part
  T4 K4  = m_K * II;

  // deviator part
  if ( epseq != 0.0 )
    K4 += 2./3. * m_sig0/std::pow(m_eps0,m_n)
       * (2./3.*(m_n-1.)*std::pow(epseq,m_n-3.)*epsd.dyadic(epsd) + std::pow(epseq,m_n-1.)*I4d);
  else
    K4 += I4d;

  return std::make_tuple(K4,sig);
}

// =================================================================================================

}}}}} // namespace ...

#endif
