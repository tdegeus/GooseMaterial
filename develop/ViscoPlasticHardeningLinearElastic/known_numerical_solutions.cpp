#include <GooseSolid/ViscoPlasticHardeningLinearElastic.h>

int main()
{

  double K      = 0.833333;
  double G      = 0.384615;
  double gamma0 = 0.01;
  double n      = 10.0;
  double sigy0  = 0.02;
  double H      = 0.08;
  double m      = 0.2;

  double dt     = 1.e-3;

  GooseSolid::ViscoPlasticHardeningLinearElastic mat(K,G,gamma0,n,sigy0,H,m);

  cppmat::tensor2s<double> eps(3,0.),sig;
  cppmat::tensor4 <double> K4(3);

  eps(0,1) = .1;
  eps(0,2) = .02;
  eps(1,2) = .045;
  eps(0,0) = .01;
  eps(1,1) = .03;
  eps(2,2) = .04;

  std::tie(K4,sig) = mat.tangent_stress(eps,dt);

  if ( std::abs( sig(0,0    ) - (  5.3846329960412831e-2 ) ) > 1.e-7 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - (  7.6921861574501663e-2 ) ) > 1.e-7 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - (  1.5384371741787811e-2 ) ) > 1.e-7 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - (  7.6921861574501663e-2 ) ) > 1.e-7 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - (  6.9230701702200650e-2 ) ) > 1.e-7 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - (  3.4614838568194528e-2 ) ) > 1.e-7 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - (  1.5384371741787811e-2 ) ) > 1.e-7 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - (  3.4614838568194528e-2 ) ) > 1.e-7 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - (  7.6922887573094553e-2 ) ) > 1.e-7 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - (  1.3461454680298226    ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - ( -2.9839651420226525e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - ( -5.9679300617225793e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - ( -2.9839651420226525e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - (  0.57692680406775565   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - ( -1.3427843472586022e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - ( -5.9679300617225793e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - ( -1.3427843472586022e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - (  0.57692677422810534   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - ( -2.9839651420226525e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - (  0.38461109252053505   ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - (  0.38461109252053505   ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - (  5.9679302840453008e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - (  2.3871721136181209e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - ( -5.9679300617225793e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - (  0.38460937375654386   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - (  1.1935860123445153e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - (  0.38460937375654386   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - (  4.7743440493780617e-8 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - ( -2.9839651420226525e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - (  0.38461109252053505   ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - (  0.38461109252053505   ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - (  5.9679302840453008e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - (  2.3871721136181209e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - (  0.57692680406775565   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - (  5.9679302840453008e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - (  1.1935860123445153e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - (  5.9679302840453008e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - (  1.3461454202863821    ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - (  2.6855686945172029e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - (  1.1935860123445153e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - (  2.6855686945172029e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - (  0.57692682197154588   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - ( -1.3427843472586022e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - (  2.6855686945172029e-8 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - (  0.38460966469317937   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - (  0.38460966469317937   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - (  1.0742274778068812e-7 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - ( -5.9679300617225793e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - (  0.38460937375654386   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - (  3.5807581704271816e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - (  1.1935860123445153e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - (  0.38460937375654386   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - (  4.7743440493780617e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - ( -1.3427843472586022e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - (  8.0567063836873088e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - (  2.6855686945172029e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - (  0.38460966469317937   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - (  1.6113412167103221e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - (  0.38460966469317937   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - (  1.0742274778068812e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - (  0.57692677422810534   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - (  2.3871721136181209e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - (  4.7743440493780617e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - (  2.3871721136181209e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - (  0.57692682197154588   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - (  1.0742274778068812e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - (  4.7743440493780617e-8 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - (  1.0742274778068812e-7 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - (  1.3461454501260324    ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;

  return 0;
}



