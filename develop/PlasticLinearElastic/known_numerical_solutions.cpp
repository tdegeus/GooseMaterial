#include <GooseSolid/PlasticLinearElastic.h>

int main()
{

// =================================================================================================

{
  double K      = 0.833333;
  double G      = 0.384615;
  double sigy0  = 0.08;
  double H      = 0.02;

  GooseSolid::PlasticLinearElastic mat(K,G,sigy0,H);

  cppmat::tensor2s<double> eps(3,0.),sig;
  cppmat::tensor4 <double> K4;

  eps(0,1) = .1;
  eps(0,2) = .02;
  eps(1,2) = .045;
  eps(0,0) = .01;
  eps(1,1) = .03;
  eps(2,2) = .04;

  std::tie(K4,sig) = mat.stress_tangent(eps);

  if ( std::abs( sig(0,0    ) - (  5.9722673410852731e-2 ) ) > 1.e-7 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - (  4.1663799558397142e-2 ) ) > 1.e-7 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - (  8.3327596012599363e-3 ) ) > 1.e-7 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - (  4.1663799558397142e-2 ) ) > 1.e-7 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - (  6.8055433012112670e-2 ) ) > 1.e-7 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - (  1.8748710266907952e-2 ) ) > 1.e-7 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - (  8.3327596012599363e-3 ) ) > 1.e-7 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - (  1.8748710266907952e-2 ) ) > 1.e-7 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - (  7.2221812812742633e-2 ) ) > 1.e-7 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - (  1.1066640713584914    ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - (  2.6565622992390101e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - (  5.3131244005487077e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - (  2.6565622992390101e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - (  0.69533920638345881   ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - (  1.1954530643469513e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - (  5.3131244005487077e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - (  1.1954530643469513e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - (  0.69799576858373324   ) ) > 1.e-7 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - (  2.6565622992390101e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - (  4.8925250795570585e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - (  4.8925250795570585e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - ( -5.3131245984780169e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - ( -2.1252498393912071e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - (  5.3131244005487077e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - (  0.20194324540713232   ) ) > 1.e-7 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - ( -1.0626248801097408e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - (  0.20194324540713232   ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - ( -4.2504995204389641e-3 ) ) > 1.e-7 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - (  2.6565622992390101e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - (  4.8925250795570585e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - (  4.8925250795570585e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - ( -5.3131245984780169e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - ( -2.1252498393912071e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - (  0.69533920638345881   ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - ( -5.3131245984780169e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - ( -1.0626248801097408e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - ( -5.3131245984780169e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - (  1.1109145708789305    ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - ( -2.3909061286939014e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - ( -1.0626248801097408e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - ( -2.3909061286939014e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - (  0.69374526906329426   ) ) > 1.e-7 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - (  1.1954530643469513e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - ( -2.3909061286939014e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - (  0.17604175994638868   ) ) > 1.e-7 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - (  0.17604175994638868   ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - ( -9.5636245147756056e-3 ) ) > 1.e-7 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - (  5.3131244005487077e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - (  0.20194324540713232   ) ) > 1.e-7 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - ( -3.1878747590868110e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - ( -1.0626248801097408e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - (  0.20194324540713232   ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - ( -4.2504995204389641e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - (  1.1954530643469513e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - ( -7.1727186532862947e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - ( -2.3909061286939014e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - (  0.17604175994638868   ) ) > 1.e-7 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - ( -1.4345436772163411e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - (  0.17604175994638868   ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - ( -9.5636245147756056e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - (  0.69799576858373324   ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - ( -2.1252498393912071e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - ( -4.2504995204389641e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - ( -2.1252498393912071e-2 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - (  0.69374526906329426   ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - ( -9.5636245147756056e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - ( -4.2504995204389641e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - ( -9.5636245147756056e-3 ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - (  1.1082580086786562    ) ) > 1.e-7 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;
}

// =================================================================================================

{
  double K      = 0.833333;
  double G      = 0.384615;
  double sigy0  = 0.08;
  double H      = 0.02;
  double m      = 0.2 ;

  GooseSolid::PlasticLinearElastic mat(K,G,sigy0,H,m);

  cppmat::tensor2s<double> eps(3,0.),sig;
  cppmat::tensor4 <double> K4;

  eps(0,1) = .1;
  eps(0,2) = .032;
  eps(1,2) = .08;
  eps(0,0) = .021;
  eps(1,1) = .017;
  eps(2,2) = .046;

  std::tie(K4,sig) = mat.stress_tangent(eps);

  if ( std::abs( sig(0,0    ) - (  6.7206685188278395e-2 ) ) > 1.e-6 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - (  3.9904122126136740e-2 ) ) > 1.e-6 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - (  1.2769319496596180e-2 ) ) > 1.e-6 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - (  3.9904122126136740e-2 ) ) > 1.e-6 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - (  6.5610520808658007e-2 ) ) > 1.e-6 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - (  3.1923296511673901e-2 ) ) > 1.e-6 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - (  1.2769319496596180e-2 ) ) > 1.e-6 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - (  3.1923296511673901e-2 ) ) > 1.e-6 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - (  7.7182715719812583e-2 ) ) > 1.e-6 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - (  1.0988361622348761    ) ) > 1.e-6 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - (  7.4904293169608619e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - (  2.3969374595587405e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - (  7.4904293169608619e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - (  0.69949532984188123   ) ) > 1.e-6 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - (  5.9923432303365031e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - (  2.3969374595587405e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - (  5.9923432303365031e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - (  0.70166755424892624   ) ) > 1.e-6 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - (  7.4904293169608619e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - (  9.2514481012145489e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - (  9.2514481012145489e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - (  1.1770673027440649e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - ( -1.9261102344401520e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - (  2.3969374595587405e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - (  0.18856317957475799   ) ) > 1.e-6 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - (  3.7666154915586931e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - (  0.18856317957475799   ) ) > 1.e-6 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - ( -6.1635529511174366e-3 ) ) > 1.e-6 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - (  7.4904293169608619e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - (  9.2514481012145489e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - (  9.2514481012145489e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - (  1.1770673027440649e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - ( -1.9261102344401520e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - (  0.69949532984188123   ) ) > 1.e-6 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - (  1.1770673027440649e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - (  3.7666154915586931e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - (  1.1770673027440649e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - (  1.0980657183850608    ) ) > 1.e-6 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - (  9.4165380711591310e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - (  3.7666154915586931e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - (  9.4165380711591310e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - (  0.70243799809874141   ) ) > 1.e-6 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - (  5.9923432303365031e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - (  9.4165380711591310e-3 ) ) > 1.e-6 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - (  0.13103669170695686   ) ) > 1.e-6 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - (  0.13103669170695686   ) ) > 1.e-6 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - ( -1.5408881301495642e-2 ) ) > 1.e-6 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - (  2.3969374595587405e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - (  0.18856317957475799   ) ) > 1.e-6 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - ( -3.4241961642704710e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - (  3.7666154915586931e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - (  0.18856317957475799   ) ) > 1.e-6 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - ( -6.1635529511174366e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - (  5.9923432303365031e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - ( -8.5604898127328577e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - (  9.4165380711591310e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - (  0.13103669170695686   ) ) > 1.e-6 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - ( -2.7393568293673801e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - (  0.13103669170695686   ) ) > 1.e-6 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - ( -1.5408881301495642e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - (  0.70166755424892624   ) ) > 1.e-6 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - ( -1.9261102344401520e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - ( -6.1635529511174366e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - ( -1.9261102344401520e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - (  0.70243799809874141   ) ) > 1.e-6 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - ( -1.5408881301495642e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - ( -6.1635529511174366e-3 ) ) > 1.e-6 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - ( -1.5408881301495642e-2 ) ) > 1.e-6 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - (  1.0958934939780161    ) ) > 1.e-6 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;
}

// =================================================================================================
  return 0;
}



