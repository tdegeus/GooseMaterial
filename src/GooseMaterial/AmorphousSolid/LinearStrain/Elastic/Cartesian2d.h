/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Linear elastic elastic material.

Suggested references
--------------------

*   The code + comments below.
*   docs/AmorphousSolid/LinearStrain/ElastoPlastic/readme.pdf

================================================================================================= */

#ifndef GOOSEMATERIAL_AMORPHOUSSOLID_LINEARSTRAIN_ELASTIC_CARTESIAN2D_H
#define GOOSEMATERIAL_AMORPHOUSSOLID_LINEARSTRAIN_ELASTIC_CARTESIAN2D_H

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <tuple>
#include <math.h>
#include <cppmat/tensor2.h>

// -------------------------------------------------------------------------------------------------

namespace GooseMaterial {
namespace AmorphousSolid {
namespace LinearStrain {
namespace Elastic {
namespace Cartesian2d {

// -------------------------------------------------------------------------------------------------

namespace cm   = cppmat::cartesian2d;
using     T2s  = cm::tensor2s<double>;
using     T2d  = cm::tensor2d<double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:
  double m_K; // bulk  modulus
  double m_G; // shear modulus

public:
  Material(){};
  Material(double K, double G);

  // compute stress or the energy at "Eps"
  T2s    stress(const T2s &Eps);
  double energy(const T2s &Eps);
};

// ===================================== IMPLEMENTATION : CORE =====================================

inline Material::Material(double K, double G)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;
}

// -------------------------------------------------------------------------------------------------

inline T2s Material::stress(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;

  // return stress tensor
  return ( m_K * epsm ) * I + m_G * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline double Material::energy(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/2.;
  T2s    Epsd = Eps - epsm*I;
  double epsd = std::sqrt(.5*Epsd.ddot(Epsd));

  // hydrostatic and deviatoric part of the energy
  double U = m_K * std::pow(epsm,2.);
  double V = m_G * std::pow(epsd,2.);

  // return total strain energy
  return U + V;
}

// =================================================================================================

}}}}} // namespace ...

#endif
