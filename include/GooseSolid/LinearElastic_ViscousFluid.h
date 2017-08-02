
#include <cppmat/tensor.h>

using T2 = cppmat::tensor2<double>;

namespace GooseSolid {

// ============================================ OVERVIEW ===========================================

class LinearElastic_ViscousFluid
{
private:
  bool   m_elas      ;  // elastic/viscous switch: true <-> elastic
  double m_K         ;  // bulk  modulus
  double m_G         ;  // shear modulus
  double m_sigy      ;  // yield stress
  double m_sigy_next ;  // 'next' yield stress applied directly after yielding
  double m_Tdamp     ;  // characteristic time scale damping
  double m_Tfluid    ;  // time duration that the material remains viscous after yielding
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
  LinearElastic_ViscousFluid(double K, double G, double sigy, double Tdamp, double Tfluid);

  // constitutive response: 'sig ( epsdot, dt )' as a function of the history
  T2   stress(const T2 &epsdot, const double dt);

  // update history
  void increment();

  // set the 'next' yield stress
  void setNextSigy(double sigy);
};

// ========================================= IMPLEMENTATION ========================================

LinearElastic_ViscousFluid::LinearElastic_ViscousFluid(
  double K, double G, double sigy , double Tdamp  , double Tfluid    ) :
  m_K(K)  , m_G(G)  , m_sigy(sigy), m_Tdamp(Tdamp), m_Tfluid(Tfluid)
{
  // set yield stress constant
  m_sigy_next = m_sigy;

  // resize tensors to "3" dimensions
  m_sigd  .resize(3);
  m_sigd_n.resize(3);

  // zero-initialize history
  m_sigd_n.zeros();
  m_sigm_n = 0.0;
  m_T_n    = 0.0;

  // initialize as elastic
  m_elas = true;
};

// =================================================================================================

void LinearElastic_ViscousFluid::increment()
{
  m_sigd_n = m_sigd;
  m_sigm_n = m_sigm;
  m_T      = m_T_n ;
}

// =================================================================================================

void LinearElastic_ViscousFluid::setNextSigy(double next)
{
  m_sigy_next = next;
}

// =================================================================================================

T2 LinearElastic_ViscousFluid::stress(const T2 &epsdot, double dt)
{
  double epsdotm,sigeq;
  T2 I,epsdotd;

  // set time
  m_T = m_T_n+dt;

  // decompose the strain in a hydrostatic and a deviatoric part
  I       = cppmat::identity2(3);
  epsdotm = epsdot.trace()/3.;
  epsdotd = epsdot - epsdotm*I;

  // compute the elastic stress (trial state)
  m_sigd  = m_sigd_n + dt*2.0*m_G*epsdotd;
  m_sigm  = m_sigm_n + dt*3.0*m_K*epsdotm;

  // if elastic -> determine whether to yield, based on the trial state
  if ( m_elas ) {
    sigeq = std::pow( 1.5 * m_sigd.ddot(m_sigd) , 0.5 );
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

};
