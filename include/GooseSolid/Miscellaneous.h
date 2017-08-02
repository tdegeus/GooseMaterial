
#include <string>
#include <tuple>
#include <cppmat/tensor.h>

using T2 = cppmat::tensor2<double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

// Convert a pair of elastic parameters to another pair: "E,nu", "K,G", "lambda,mu"
// E.g.: std::tie(K,G) = GooseSolid::ConvertElasticParameters("E,nu",1.,.3,"K,G");
std::tuple<double,double> ConvertElasticParameters (
  std::string in, double ipar1, double ipar2, std::string out );

// Von-Mises equivalent stress
double VonMisesStress(const T2 &stress);

// ========================================= IMPLEMENTATION ========================================

std::tuple<double,double> ConvertElasticParameters (
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
    throw std::runtime_error("GooseSolid::ConvertElasticParameters -> Unknown input pair");
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

  throw std::runtime_error("GooseSolid::ConvertElasticParameters -> Unknown output pair");
}

// =================================================================================================

double VonMisesStress(const T2 &S)
{
  // S        = sigma - trace(sigma)/3
  // sigma_eq = S_{ij} * S_{ji}
  if ( S.ndim() == 3 )
  {
    return std::pow(
      .5*( std::pow(S(0,0)-S(1,1),2.) + std::pow(S(1,1)-S(2,2),2.) + std::pow(S(2,2)-S(0,0),2.) ) +
      3.*( std::pow(S(0,1),2.) + std::pow(S(0,2),2.) + std::pow(S(1,2),2.) )
    ,0.5);
  }

  throw std::runtime_error("Unknown stress definition");
}

};
