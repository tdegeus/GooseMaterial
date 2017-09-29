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

#include <tuple>
#include <cppmat/tensor3.h>

#warning "GooseMaterial/Metal/LinearStrain/Elastic/Cartesian3d.h : first usage, careful check then remove this message"

namespace GooseMaterial {
namespace Metal {
namespace LinearStrain {
namespace Elastic {

using T2  = cppmat::tensor3_2 <double>;
using T2s = cppmat::tensor3_2s<double>;
using T2d = cppmat::tensor3_2d<double>;
using T4  = cppmat::tensor3_4 <double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:

  double m_K; // material parameter : bulk  modulus
  double m_G; // material parameter : shear modulus

  // compute the stress (and optionally the tangent)
  std::tuple<T4,T2s> compute(const T2s &eps, bool tangent);

public:

  // constructor / destructor
 ~Material(){};
  Material(){};
  Material(double K, double G);

  // compute stress(+tangent) at "eps"
  T2s                stress        (const T2s &eps);
  std::tuple<T4,T2s> tangent_stress(const T2s &eps);

};

// =================================================================================================

std::tuple<double,double> ConvertParameters (
  std::string in, double ipar1, double ipar2, std::string out );

// ========================================= IMPLEMENTATION ========================================

Material::Material( double K, double G ) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

T2s  Material::stress(const T2s &eps)
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

std::tuple<T4,T2s> Material::tangent_stress(const T2s &eps)
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

// =================================================================================================

std::tuple<double,double> ConvertParameters (
  std::string in, double ipar1, double ipar2, std::string out )
{
  double E,nu,K,G;

  // convert input to "K" and "G"
  if      ( in == "E,nu" ) {
    E  = ipar1;
    nu = ipar2;
    K  = E / ( 3.*(1.-2.*nu) );
    G  = E / ( 2.*(1.+   nu) );
  }
  else if ( in == "lambda,mu" ) {
    K  = ipar1 + 2./3.*ipar2;
    G  = ipar2;
  }
  else if ( in == "K,G" ) {
    K  = ipar1;
    G  = ipar2;
  }
  else {
    throw std::runtime_error("ConvertParameters -> Unknown input pair");
  }

  // return requested pair
  if ( out == "K,G" )
    return std::make_tuple(K,G);

  if ( out == "lambda,mu" )
    return std::make_tuple(K-2./3.*G,G);

  if ( out == "E,nu" ) {
    E  = ( 9.*K   *G ) / (     3.*K+G  );
    nu = ( 3.*K-2.*G ) / ( 2.*(3.*K+G) );
    return std::make_tuple(E,nu);
  }

  throw std::runtime_error("ConvertParameters -> Unknown output pair");
}

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
