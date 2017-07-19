
#include <cppmat/tensor.h>

using T2 = cppmat::tensor2<double>;
using T4 = cppmat::tensor4<double>;

namespace GooseSolid {

// =================================================================================================

class ViscoPlasticLinearElastic
{
private:

public:
 ~ViscoPlasticLinearElastic(){};
  ViscoPlasticLinearElastic(){};
  ViscoPlasticLinearElastic(
    double K, double mu
  );

  T4     tangent(bool stress_only=false);

  T2     stress();

  T2     getStress();

  void   increment();

};

// =================================================================================================

ViscoPlasticLinearElastic::ViscoPlasticLinearElastic(
  double K, double mu ) :
  m_K(K)  , m_mu(mu)
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

void ViscoPlasticLinearElastic::setNext_sigy(double next)
{
  m_sigy_next = next;
}

// =================================================================================================

double ViscoPlasticLinearElastic::getEquivStress()
{
  return pow(0.5*m_sigd.ddot(m_sigd),0.5);
}

// =================================================================================================

T2 ViscoPlasticLinearElastic::stress(T2 epsdot, double dt)
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

void ViscoPlasticLinearElastic::increment()
{
  m_sigd_n = m_sigd;
  m_sigm_n = m_sigm;
  m_T      = m_T_n ;
}

// =================================================================================================

};
