/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Linear elasticity: a linear relationship between the linear strain and the Cauchy stress.

Suggested references
--------------------

*   The code + comments below.
*   docs/Metal/LinearStrain/Elastic/readme.pdf
*   Former internal code: GooseFEM / mat1001

================================================================================================= */

#ifndef GOOSEMATERIAL_METAL_LINEARSTRAIN_ELASTIC_CARTESIAN3D_H
#define GOOSEMATERIAL_METAL_LINEARSTRAIN_ELASTIC_CARTESIAN3D_H

#include <tuple>
#include <cppmat/tensor3.h>

#include "../../../Macros.h"

// -------------------------------------------------------------------------------------------------

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace Elastic {
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
  double m_K; // material parameter : bulk  modulus
  double m_G; // material parameter : shear modulus

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

public:
  Material(){};
  Material(double K, double G);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

Material::Material( double K, double G ) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

T2s  Material::stress(const T2s &eps)
{
  // second order identity tensor
  T2d    I    = cm::identity2();

  // decompose strain: hydrostatic part, deviatoric part
  double epsm = eps.trace() / 3.;
  T2s    epsd = eps - epsm * I;

  // constitutive response
  double sigm = ( 3. * m_K ) * epsm;
  T2s    sigd = ( 2. * m_G ) * epsd;

  // combine volumetric and deviatoric stress
  return sigm * I + sigd ;
}

// -------------------------------------------------------------------------------------------------

std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps)
{
  // stress
  // ------

  // second order identity tensor
  T2d    I    = cm::identity2();

  // decompose strain: hydrostatic part, deviatoric part
  double epsm = eps.trace() / 3.;
  T2s    epsd = eps - epsm * I;

  // constitutive response
  double sigm = ( 3. * m_K ) * epsm;
  T2s    sigd = ( 2. * m_G ) * epsd;

  // combine volumetric and deviatoric stress
  T2s    sig  = sigm * I + sigd;

  // tangent
  // -------

  // unit tensors: II = dyadic(I,I) and deviatoric unit tensor I4d (A_d = I4d : A)
  T4 I4d = cm::identity4d();
  T4 II  = cm::identityII();

  // initialize tangent as the elasticity tensor
  T4 K4  = m_K * II + ( 2. * m_G ) * I4d;

  return std::make_tuple(K4,sig);
}

// =================================================================================================

}}}}} // namespace GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d {

#endif
