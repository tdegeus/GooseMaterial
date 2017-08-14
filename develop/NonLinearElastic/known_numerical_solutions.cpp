#include <GooseSolid/NonLinearElastic.h>

int main()
{

// =================================================================================================

{
  double K     =  0.833333;
  double sig0  =  0.5;
  double eps0  =  0.1;
  double n     = 10.1;

  GooseSolid::NonLinearElastic mat(K,sig0,eps0,n);

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

  if ( std::abs( sig(0,0    ) - ( -1.0820895138976112 ) ) > 1.e-6 ) std::cout << "sig(0,0    ) failed" << std::endl;
  if ( std::abs( sig(0,1    ) - (  16.458420261654521 ) ) > 1.e-6 ) std::cout << "sig(0,1    ) failed" << std::endl;
  if ( std::abs( sig(0,2    ) - (  5.2666946554041445 ) ) > 1.e-6 ) std::cout << "sig(0,2    ) failed" << std::endl;
  if ( std::abs( sig(1,0    ) - (  16.458420261654521 ) ) > 1.e-6 ) std::cout << "sig(1,0    ) failed" << std::endl;
  if ( std::abs( sig(1,1    ) - ( -1.7404261159016579 ) ) > 1.e-6 ) std::cout << "sig(1,1    ) failed" << std::endl;
  if ( std::abs( sig(1,2    ) - (  13.166735718824476 ) ) > 1.e-6 ) std::cout << "sig(1,2    ) failed" << std::endl;
  if ( std::abs( sig(2,0    ) - (  5.2666946554041445 ) ) > 1.e-6 ) std::cout << "sig(2,0    ) failed" << std::endl;
  if ( std::abs( sig(2,1    ) - (  13.166735718824476 ) ) > 1.e-6 ) std::cout << "sig(2,1    ) failed" << std::endl;
  if ( std::abs( sig(2,2    ) - (  3.0325155515160187 ) ) > 1.e-6 ) std::cout << "sig(2,2    ) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,0) - (  112.63264573142889 ) ) > 1.e-4 ) std::cout << "K4 (0,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,1) - ( -29.664464019677503 ) ) > 1.e-4 ) std::cout << "K4 (0,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,0,2) - ( -9.4926287957212701 ) ) > 1.e-4 ) std::cout << "K4 (0,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,0) - ( -29.664464019677503 ) ) > 1.e-4 ) std::cout << "K4 (0,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,1) - ( -50.764976247563219 ) ) > 1.e-4 ) std::cout << "K4 (0,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,1,2) - ( -23.731570331672096 ) ) > 1.e-4 ) std::cout << "K4 (0,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,0) - ( -9.4926287957212701 ) ) > 1.e-4 ) std::cout << "K4 (0,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,1) - ( -23.731570331672096 ) ) > 1.e-4 ) std::cout << "K4 (0,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,0,2,2) - ( -59.367670437539985 ) ) > 1.e-4 ) std::cout << "K4 (0,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,0) - ( -29.664464019677503 ) ) > 1.e-4 ) std::cout << "K4 (0,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,1) - (  506.07013194614655 ) ) > 1.e-4 ) std::cout << "K4 (0,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,0,2) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (0,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,0) - (  506.07013194614655 ) ) > 1.e-4 ) std::cout << "K4 (0,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,1) - ( -46.615579926675409 ) ) > 1.e-4 ) std::cout << "K4 (0,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,1,2) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (0,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,0) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (0,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,1) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (0,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,1,2,2) - (  76.280043946352919 ) ) > 1.e-4 ) std::cout << "K4 (0,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,0) - ( -9.4926287957212701 ) ) > 1.e-4 ) std::cout << "K4 (0,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,1) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (0,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,0,2) - (  125.68697337393439 ) ) > 1.e-4 ) std::cout << "K4 (0,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,0) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (0,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,1) - ( -14.916986062774512 ) ) > 1.e-4 ) std::cout << "K4 (0,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,1,2) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (0,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,0) - (  125.68697337393439 ) ) > 1.e-4 ) std::cout << "K4 (0,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,1) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (0,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (0,2,2,2) - (  24.409614858495786 ) ) > 1.e-4 ) std::cout << "K4 (0,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,0) - ( -29.664464019677503 ) ) > 1.e-4 ) std::cout << "K4 (1,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,1) - (  506.07013194614655 ) ) > 1.e-4 ) std::cout << "K4 (1,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,0,2) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (1,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,0) - (  506.07013194614655 ) ) > 1.e-4 ) std::cout << "K4 (1,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,1) - ( -46.615579926675409 ) ) > 1.e-4 ) std::cout << "K4 (1,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,1,2) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (1,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,0) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (1,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,1) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (1,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,0,2,2) - (  76.280043946352919 ) ) > 1.e-4 ) std::cout << "K4 (1,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,0) - ( -50.764976247563219 ) ) > 1.e-4 ) std::cout << "K4 (1,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,1) - ( -46.615579926675409 ) ) > 1.e-4 ) std::cout << "K4 (1,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,0,2) - ( -14.916986062774512 ) ) > 1.e-4 ) std::cout << "K4 (1,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,0) - ( -46.615579926675409 ) ) > 1.e-4 ) std::cout << "K4 (1,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,1) - (  115.68384652312099 ) ) > 1.e-4 ) std::cout << "K4 (1,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,1,2) - ( -37.292462552087805 ) ) > 1.e-4 ) std::cout << "K4 (1,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,0) - ( -14.916986062774512 ) ) > 1.e-4 ) std::cout << "K4 (1,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,1) - ( -37.292462552087805 ) ) > 1.e-4 ) std::cout << "K4 (1,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,1,2,2) - ( -62.418871229232067 ) ) > 1.e-4 ) std::cout << "K4 (1,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,0) - ( -23.731570331672096 ) ) > 1.e-4 ) std::cout << "K4 (1,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,1) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (1,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,0,2) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (1,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,0) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (1,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,1) - ( -37.292462552087805 ) ) > 1.e-4 ) std::cout << "K4 (1,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,1,2) - (  353.51002026775217 ) ) > 1.e-4 ) std::cout << "K4 (1,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,0) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (1,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,1) - (  353.51002026775217 ) ) > 1.e-4 ) std::cout << "K4 (1,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (1,2,2,2) - (  61.024032883759922 ) ) > 1.e-4 ) std::cout << "K4 (1,2,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,0) - ( -9.4926287957212701 ) ) > 1.e-4 ) std::cout << "K4 (2,0,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,1) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (2,0,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,0,2) - (  125.68697337393439 ) ) > 1.e-4 ) std::cout << "K4 (2,0,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,0) - (  135.60897461686824 ) ) > 1.e-4 ) std::cout << "K4 (2,0,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,1) - ( -14.916986062774512 ) ) > 1.e-4 ) std::cout << "K4 (2,0,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,1,2) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (2,0,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,0) - (  125.68697337393439 ) ) > 1.e-4 ) std::cout << "K4 (2,0,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,1) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (2,0,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,0,2,2) - (  24.409614858495786 ) ) > 1.e-4 ) std::cout << "K4 (2,0,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,0) - ( -23.731570331672096 ) ) > 1.e-4 ) std::cout << "K4 (2,1,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,1) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (2,1,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,0,2) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (2,1,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,0) - (  339.02241286172807 ) ) > 1.e-4 ) std::cout << "K4 (2,1,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,1) - ( -37.292462552087805 ) ) > 1.e-4 ) std::cout << "K4 (2,1,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,1,2) - (  353.51002026775217 ) ) > 1.e-4 ) std::cout << "K4 (2,1,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,0) - (  108.48717565203226 ) ) > 1.e-4 ) std::cout << "K4 (2,1,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,1) - (  353.51002026775217 ) ) > 1.e-4 ) std::cout << "K4 (2,1,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,1,2,2) - (  61.024032883759922 ) ) > 1.e-4 ) std::cout << "K4 (2,1,2,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,0) - ( -59.367670437539985 ) ) > 1.e-4 ) std::cout << "K4 (2,2,0,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,1) - (  76.280043946352919 ) ) > 1.e-4 ) std::cout << "K4 (2,2,0,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,0,2) - (  24.409614858495786 ) ) > 1.e-4 ) std::cout << "K4 (2,2,0,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,0) - (  76.280043946352919 ) ) > 1.e-4 ) std::cout << "K4 (2,2,1,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,1) - ( -62.418871229232067 ) ) > 1.e-4 ) std::cout << "K4 (2,2,1,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,1,2) - (  61.024032883759922 ) ) > 1.e-4 ) std::cout << "K4 (2,2,1,2) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,0) - (  24.409614858495786 ) ) > 1.e-4 ) std::cout << "K4 (2,2,2,0) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,1) - (  61.024032883759922 ) ) > 1.e-4 ) std::cout << "K4 (2,2,2,1) failed" << std::endl;
  if ( std::abs( K4 (2,2,2,2) - (  124.28654071309776 ) ) > 1.e-4 ) std::cout << "K4 (2,2,2,2) failed" << std::endl;
}




























































































// =================================================================================================

  return 0;
}



