
#include <cppmat/tensor.h>

using T2 = cppmat::tensor2<double>;

namespace GooseSolid {

// =================================================================================================

class LinearElastic_ViscousFluid
{
private:
  bool   m_elas      ;  // elastic/viscous switch: true <-> elastic
  double m_K         ;  // bulk  modulus
  double m_mu        ;  // shear modulus
  double m_sigy      ;  // yield stress
  double m_sigy_next ;  // 'next' yield stress applied directly after yielding
  double m_Tdamp     ;  // characteristic time scale damping
  double m_Tfluid    ;  // time duration that the material remains viscous after yielding
  size_t m_nd        ;  // number of dimensions
  double m_T         ;  // current time
  double m_T_n       ;  // time of the last increment
  double m_Tyield_n  ;  // time at which the material has yielded last
  T2     m_sigd      ;  // deviatoric  stress
  T2     m_sigd_n    ;  // deviatoric  stress at previous time-step
  double m_sigm      ;  // hydrostatic stress
  double m_sigm_n    ;  // hydrostatic stress at previous time-step

public:
 ~LinearElastic_ViscousFluid(){};
  LinearElastic_ViscousFluid(){};
  LinearElastic_ViscousFluid(
    double K, double mu, double sigy, double Tdamp, double Tfluid, size_t nd=3
  );

  // constitutive response: 'sigma_{n+1} ( epsdot_{n+1}, dt )'
  T2     stress(T2 epsdot, double dt);
  // return (equivalent) stress (to get it at 'n+1', call after calling "stress" for that increment)
  T2     getStress();
  double getEquivStress();
  // set the 'next' yield stress
  void   setNext_sigy(double sigy);
  // store history: 'sigd', 'sigm' -> 'sigd_n', 'sigm_n'
  void   increment();

};

// =================================================================================================

LinearElastic_ViscousFluid::LinearElastic_ViscousFluid(
  double K, double mu, double sigy , double Tdamp  , double Tfluid   , size_t nd ) :
  m_K(K)  , m_mu(mu) , m_sigy(sigy), m_Tdamp(Tdamp), m_Tfluid(Tfluid), m_nd(nd)
{
  // set yield stress constant
  m_sigy_next = m_sigy;

  // resize tensors to "m_nd" dimensions
  m_sigd  .resize(m_nd);
  m_sigd_n.resize(m_nd);

  // zero-initialize history
  m_sigd_n.zeros();
  m_sigm_n = 0.0;
  m_T_n    = 0.0;

  // initialize as elastic
  m_elas = true;
};

// =================================================================================================

void LinearElastic_ViscousFluid::setNext_sigy(double next)
{
  m_sigy_next = next;
}

// =================================================================================================

T2 LinearElastic_ViscousFluid::getStress()
{
  T2 I = cppmat::identity2(m_nd);

  return m_sigd+m_sigm*I;
}

// =================================================================================================

double LinearElastic_ViscousFluid::getEquivStress()
{
  return pow(0.5*m_sigd.ddot(m_sigd),0.5);
}

// =================================================================================================

T2 LinearElastic_ViscousFluid::stress(T2 epsdot, double dt)
{
  double epsdotm;
  T2 I,epsdotd;

  // set time
  m_T = m_T_n+dt;

  // decompose the strain in a hydrostatic and a deviatoric part
  I       = cppmat::identity2(m_nd);
  epsdotm = epsdot.trace()/static_cast<double>(m_nd);
  epsdotd = epsdot-epsdotm*I;

  // compute the elastic stress (trial state)
  m_sigd  = m_sigd_n+dt*m_mu*epsdotd;
  m_sigm  = m_sigm_n+dt*m_K *epsdotm;

  // if elastic -> determine whether to yield, based on the trial state
  if ( m_elas ) {
    double sigeq = getEquivStress();
    if ( sigeq >= m_sigy ) {
      m_elas     = false;       // switch to plastic response
      m_Tyield_n = m_T;         // store current time, end yield at "m_Tyield_n + m_Tfluid"
      m_sigy     = m_sigy_next; // set next yield stress: will play a role once elastic again
    }
  }

  // if plastic -> determine whether to end yielding, based on time criterion
  if ( !m_elas )
    if ( m_T > m_Tyield_n+m_Tfluid )
      m_elas = true;

  // if plastic -> correct trial state (deviator only)
  if ( !m_elas )
    m_sigd *= m_Tdamp/(m_Tdamp+dt);

  // return total stress tensor
  return m_sigd+m_sigm*I;
};

// =================================================================================================

void LinearElastic_ViscousFluid::increment()
{
  m_sigd_n = m_sigd;
  m_sigm_n = m_sigm;
  m_T      = m_T_n ;
}

// =================================================================================================

};
