#include <GooseSolid/ViscoPlasticLinearElastic.h>

int main()
{

  double K      = 0.833333;
  double G      = 0.384615;
  double sig0   = 0.01;
  double m      = 10.0;
  double gamma0 = 1.0;
  double dt     = 1.e-3;

  GooseSolid::ViscoPlasticLinearElastic mat(K,G,sig0,gamma0,m);

  cppmat::tensor2s<double> eps(3,0.),sig;
  cppmat::tensor4 <double> K4(3);

  eps(0,1) = .1;
  eps(0,2) = .02;
  eps(1,2) = .045;
  eps(0,0) = .01;
  eps(1,1) = .03;
  eps(2,2) = .04;

  std::tie(K4,sig) = mat.tangent_stress(eps,dt);

  if ( std::abs( sig(0,0    ) - (  5.3975377048306072e-2 ) ) > 1.e-7 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - (  7.6147579018297926e-2 ) ) > 1.e-7 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - (  1.5229515236315921e-2 ) ) > 1.e-7 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - (  7.6147579018297926e-2 ) ) > 1.e-7 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - (  6.9204892284621997e-2 ) ) > 1.e-7 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - (  3.4266411409249568e-2 ) ) > 1.e-7 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - (  1.5229515236315921e-2 ) ) > 1.e-7 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - (  3.4266411409249568e-2 ) ) > 1.e-7 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - (  7.6819649902779952e-2 ) ) > 1.e-7 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - (  1.3410600290536538    ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - ( -4.5896634323029835e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - ( -9.1793265226493986e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - ( -4.5896634323029835e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - (  0.57949245695232154   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - ( -2.0653485958298278e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - ( -9.1793265226493986e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - ( -2.0653485958298278e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - (  0.57944656031970820   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - ( -4.5896634323029835e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - (  0.38349168758002172   ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - (  0.38349168758002172   ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - (  9.1793268646059615e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - (  3.6717307458423851e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - ( -9.1793265226493986e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - (  0.38084804133632477   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - (  1.8358653045298786e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - (  0.38084804133632477   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - (  7.3434612181195156e-5 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - ( -4.5896634323029835e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - (  0.38349168758002172   ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - (  0.38349168758002172   ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - (  9.1793268646059615e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - (  3.6717307458423851e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - (  0.57949245695232154   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - (  9.1793268646059615e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - (  1.8358653045298786e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - (  9.1793268646059615e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - (  1.3409865944414727    ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - (  4.1306971916596528e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - (  1.8358653045298786e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - (  4.1306971916596528e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - (  0.57951999493188944   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - ( -2.0653485958298278e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - (  4.1306971916596528e-5 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - (  0.38129553357355012   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - (  0.38129553357355012   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - (  1.6522788766638611e-4 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - ( -9.1793265226493986e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - (  0.38084804133632477   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - (  5.5075961187635780e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - (  1.8358653045298786e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - (  0.38084804133632477   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - (  7.3434612181195156e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - ( -2.0653485958298278e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - (  1.2392092036620358e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - (  4.1306971916596528e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - (  0.38129553357355012   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - (  2.4784183149957921e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - (  0.38129553357355012   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - (  1.6522788766638611e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - (  0.57944656031970820   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - (  3.6717307458423851e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - (  7.3434612181195156e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - (  3.6717307458423851e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - (  0.57951999493188944   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - (  1.6522788766638611e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - (  7.3434612181195156e-5 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - (  1.6522788766638611e-4 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - (  1.3410324910740858    ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;

  return 0;
}



