/* ========================================== DESCRIPTION ==========================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseSolid

Overview
--------

- ConvertElasticParameters    // convert pairs of elastic parameters

Description
-----------

Miscellaneous support functions.

Suggested references
--------------------

*   The code + comments below.

================================================================================================= */

#include <string>
#include <tuple>
#include <cppmat/tensor.h>

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

// Convert a pair of elastic parameters to another pair: "E,nu", "K,G", "lambda,mu". For example:
//
//    std::tie(K,G) = GooseSolid::ConvertElasticParameters("E,nu",1.,.3,"K,G");

std::tuple<double,double> ConvertElasticParameters (
  std::string in, double ipar1, double ipar2, std::string out );

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

// -------------------------------------------------------------------------------------------------

}
