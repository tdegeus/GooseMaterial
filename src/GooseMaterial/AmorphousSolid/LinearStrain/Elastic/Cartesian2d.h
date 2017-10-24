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

#include "../../../Macros.h"

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
double    ndim = 2.;

// ============================================ OVERVIEW ===========================================

class Material
{
private:
  double m_K; // bulk  modulus
  double m_G; // shear modulus

public:
  Material(){};
  Material(double K, double G);

  // compute stress at "Eps"
  T2s stress(const T2s &Eps);

  // post-process functions
  // - return (hydrostatic/deviatoric) equivalent stress/strain
  double eps_eq(const T2s &Eps);
  double eps_m (const T2s &Eps);
  double eps_d (const T2s &Eps);
  double sig_eq(const T2s &Sig);
  double sig_m (const T2s &Sig);
  double sig_d (const T2s &Sig);
  // - return the strain energy (or its hydrostatic or deviatoric component)
  double energy  (double epsm, double epsd);
  double energy_m(double epsm);
  double energy_d(double epsd);
  double energy  (const T2s &Eps);
  double energy_m(const T2s &Eps);
  double energy_d(const T2s &Eps);

};

// ===================================== IMPLEMENTATION : CORE =====================================

Material::Material(double K, double G)
{
  // copy input - elastic moduli
  m_K = K;
  m_G = G;
}

// -------------------------------------------------------------------------------------------------

T2s Material::stress(const T2s &Eps)
{
  // decompose strain: hydrostatic part, deviatoric part
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/ndim;
  T2s    Epsd = Eps - epsm*I;

  // return stress tensor
  return ( m_K * epsm ) * I + m_G * Epsd;
}

// ================================= IMPLEMENTATION : POST-PROCESS =================================

double Material::eps_eq(const T2s &Eps)
{
  return std::pow( .5*Eps.ddot(Eps) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::eps_m(const T2s &Eps)
{
  return Eps.trace()/ndim;
}

// -------------------------------------------------------------------------------------------------

double Material::eps_d(const T2s &Eps)
{
  T2d    I    = cm::identity2();
  double epsm = Eps.trace()/ndim;
  T2s    Epsd = Eps - epsm*I;

  return std::pow( .5*Epsd.ddot(Epsd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_eq(const T2s &Sig)
{
  return std::pow( .5*Sig.ddot(Sig) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_m(const T2s &Sig)
{
  return Sig.trace()/ndim;
}

// -------------------------------------------------------------------------------------------------

double Material::sig_d(const T2s &Sig)
{
  T2d    I    = cm::identity2();
  double sigm = Sig.trace()/ndim;
  T2s    Sigd = Sig - sigm*I;

  return std::pow( .5*Sigd.ddot(Sigd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(double epsm)
{
  return ndim/2. * m_K * std::pow( epsm , 2. );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_m(const T2s &Eps)
{
  return energy_m(eps_m(Eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy_d(double epsd)
{
  return m_G * std::pow( epsd , 2. );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_d(const T2s &Eps)
{
  return energy_d(eps_d(Eps));
}

// -------------------------------------------------------------------------------------------------

double Material::energy(double epsm, double epsd)
{
  return energy_m(epsm) + energy_d(epsd);
}

// -------------------------------------------------------------------------------------------------

double Material::energy(const T2s &Eps)
{
  return energy_m(eps_m(Eps)) + energy_d(eps_d(Eps));
}

// =================================================================================================

}}}}} // namespace GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian2d

#endif
