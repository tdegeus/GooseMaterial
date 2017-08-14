#include <GooseSolid/LinearElastic.h>

int main()
{

// =================================================================================================

{
  double K      = 0.833333;
  double G      = 0.384615;

  GooseSolid::LinearElastic mat(K,G);

  cppmat::tensor2s<double> eps(3,0.),sig;
  cppmat::tensor4 <double> K4;

  eps(0,1) = .1;
  eps(0,2) = .032;
  eps(1,2) = .08;
  eps(0,0) = .021;
  eps(1,1) = .017;
  eps(2,2) = .046;

  std::tie(K4,sig) = mat.tangent_stress(eps);

  std::cout << "Test 1 (no output -> completed successfully)" << std::endl;

  if ( std::abs( sig(0,0    ) - ( 6.4615363443700335e-2 ) ) > 1.e-7 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - ( 7.6923001958772552e-2 ) ) > 1.e-7 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - ( 2.4615361429176641e-2 ) ) > 1.e-7 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - ( 7.6923001958772552e-2 ) ) > 1.e-7 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - ( 6.1538444339655163e-2 ) ) > 1.e-7 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - ( 6.1538399274533973e-2 ) ) > 1.e-7 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - ( 2.4615361429176641e-2 ) ) > 1.e-7 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - ( 6.1538399274533973e-2 ) ) > 1.e-7 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - ( 8.3846113933393474e-2 ) ) > 1.e-7 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - ( 1.3461530208587646    ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - ( 1.3461530208587646    ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - ( 0.38461500406265259   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - ( 0.57692301273345947   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - ( 0.0                   ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - ( 1.3461530208587646    ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;
}

// =================================================================================================

  return 0;
}



