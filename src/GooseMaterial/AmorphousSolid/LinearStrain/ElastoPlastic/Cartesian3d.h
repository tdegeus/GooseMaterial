/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseMaterial

Description
-----------

Elastic plastic material for amorphous solids.

Suggested references
--------------------

*   The code + comments below.
*   docs/AmorphousSolid/LinearStrain/ElastoPlastic/readme.pdf

================================================================================================= */

#define _USE_MATH_DEFINES // to use "M_PI" from "math.h"

#include <tuple>
#include <math.h>
#include <cppmat/tensor3.h>

#warning "GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian3d.h : first usage, careful check then remove this message"

namespace GooseMaterial {
namespace AmorphousSolid {
namespace LinearStrain {
namespace ElastoPlastic {

using T2  = cppmat::tensor3_2 <double>;
using T2s = cppmat::tensor3_2s<double>;
using T2d = cppmat::tensor3_2d<double>;

// ============================================ OVERVIEW ===========================================

class Material
{
private:

  double              m_K;      // bulk  modulus
  double              m_G;      // shear modulus
  std::vector<double> m_epsy;   // yield strains
  bool                m_smooth; // set potential smooth or not

public:

  // constructor / destructor
 ~Material(){};
  Material(){};
  Material(double K, double G,
    const std::vector<double> &epsy={}, bool init_elastic=true, bool smooth=true);

  // compute stress at "eps"
  T2s stress(const T2s &eps);

  // post-process functions
  // - find the index of the current yield strain (below "eps_eq")
  size_t find(double epseq);
  // - return yield strain "i"
  double eps_y(size_t i);
  // - return equivalent strain
  double eps_eq(const T2s &eps);
  // - return hydrostatic stress
  double sig_m(const T2s &sig);
  // - return equivalent stress
  double sig_eq(const T2s &sig);
  // - return the strain energy related to the equivalent strains
  double energy_eq(const T2s &eps);

};

// ========================================= IMPLEMENTATION ========================================

Material::Material(
  double K, double G, const std::vector<double> &epsy, bool init_elastic, bool smooth )
{
  m_K      = K;
  m_G      = G;
  m_smooth = smooth;

  // copy input
  std::vector<double> vec = epsy;
  // initialize
  if ( vec.size() == 0 ) vec = { 1000. };
  // sort input
  std::sort( vec.begin() , vec.end() );

  // check to add item to have an initial elastic response
  if ( init_elastic )
    if ( vec[0] == -vec[1] )
      init_elastic = false;

  // copy input
  // - counters
  size_t N = vec.size();
  size_t i = 0;
  // - add yield stress to have an initial elastic response
  if ( init_elastic ) ++N;
  // - allocate
  m_epsy.resize(N);
  // - add yield stress to have an initial elastic response
  if ( init_elastic ) { m_epsy[i] = -vec[0]; ++i; }
  // - copy the rest
  for ( auto &j : vec ) {
    m_epsy[i] = j; ++i;
  }

  // check the number of yield strains
  if ( m_epsy.size() < 2 )
    throw std::runtime_error("Specify at least two yield strains 'eps_y'");
}

// -------------------------------------------------------------------------------------------------

size_t Material::find(double epseq)
{
  // check extremes
  if ( epseq < m_epsy.front() or epseq >= m_epsy.back() )
    throw std::runtime_error("Insufficient 'eps_y'");

  // set initial search bounds and index
  size_t n = m_epsy.size()-1;  // upper-bound
  size_t z = 1;                // lower-bound
  size_t l = 0;                // left-bound
  size_t r = n;                // right-bound
  size_t i = r / 2;            // average

  // loop until found
  while ( true )
  {
    // check if found, unroll once to speed-up
    if ( epseq >= m_epsy[i-1] and epseq < m_epsy[i  ] ) return i-1;
    if ( epseq >= m_epsy[i  ] and epseq < m_epsy[i+1] ) return i;
    if ( epseq >= m_epsy[i+1] and epseq < m_epsy[i+2] ) return i+1;

    // correct the left- and right-bound
    if ( epseq >= m_epsy[i] ) l = i;
    else                      r = i;

    // set new search index
    i = ( r + l ) / 2;
    i = std::max(i,z);
  }
}

// -------------------------------------------------------------------------------------------------

double Material::eps_eq(const T2s &eps)
{
  double epsm;
  T2s epsd;
  T2d I;

  // second order identity tensor
  I    = cppmat::identity3_2();

  // decompose strain: hydrostatic part, deviatoric part
  epsm = eps.trace() / 3.;
  epsd = eps - epsm * I;

  // equivalent strain
  return std::pow( 2./3. * epsd.ddot(epsd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::sig_m(const T2s &sig)
{
  return sig.trace() / 3.;
}

// -------------------------------------------------------------------------------------------------

double Material::sig_eq(const T2s &sig)
{
  double sigm;
  T2s sigd;
  T2d I;

  // second order identity tensor
  I    = cppmat::identity3_2();

  // decompose stress: hydrostatic part, deviatoric part
  sigm = sig.trace() / 3.;
  sigd = sig - sigm * I;

  // equivalent stress
  return std::pow( 3./2. * sigd.ddot(sigd) , 0.5 );
}

// -------------------------------------------------------------------------------------------------

double Material::energy_eq(const T2s &eps)
{
  double epsm,epseq;
  T2s epsd;
  T2d I;

  // second order identity tensor
  I     = cppmat::identity3_2();

  // decompose strain: hydrostatic part, deviatoric part
  epsm  = eps.trace() / 3.;
  epsd  = eps - epsm * I;

  // equivalent strain
  epseq = std::pow( 2./3. * epsd.ddot(epsd) , 0.5 );

  // read current yield strains
  size_t i       = find(epseq);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // energy
  if ( ! m_smooth )
    return 1.5 * m_G * ( std::pow( epseq-eps_min , 2. ) - std::pow( deps_y , 2. ) );
  else
    return -3. * m_G * std::pow( deps_y/M_PI , 2. ) * ( 1. + cos( M_PI/deps_y * (epseq-eps_min) ) );
}

// -------------------------------------------------------------------------------------------------

double Material::eps_y(size_t i)
{
  return m_epsy[i];
}

// -------------------------------------------------------------------------------------------------

T2s  Material::stress(const T2s &eps)
{
  double epsm,sigm,epseq;
  T2s epsd,sigd,sig;
  T2d I;

  // second order identity tensor
  I     = cppmat::identity3_2();

  // decompose strain: hydrostatic part, deviatoric part
  epsm  = eps.trace() / 3.;
  epsd  = eps - epsm * I;

  // equivalent strain
  epseq = std::pow( 2./3. * epsd.ddot(epsd) , 0.5 );

  // constitutive response - hydrostatic part
  sigm  = 3. * m_K * epsm;

  // equivalent strain zero -> zero deviatoric stress
  if ( epseq <= 0. ) return sigm * I;

  // read current yield strains
  size_t i       = find(epseq);
  double eps_min = ( m_epsy[i+1] + m_epsy[i] ) / 2.;
  double deps_y  = ( m_epsy[i+1] - m_epsy[i] ) / 2.;

  // constitutive response - deviatoric part
  if ( ! m_smooth )
    sigd = 2.*m_G * ( 1. - eps_min/epseq ) * epsd;
  else
    sigd = 2.*m_G * ( deps_y/M_PI ) * sin ( M_PI/deps_y * ( epseq - eps_min ) ) * epsd / epseq;

  // return full strain tensor
  return sigm * I + sigd;
}

// =================================================================================================

} // namespace ...
} // namespace ...
} // namespace ...
} // namespace ...
