
#include <cppmat/tensor.h>
#include <cppmat/pybind11_tensor.h>

#include "include/GooseSolid/LinearElastic_ViscousFluid.h"
#include "include/GooseSolid/ViscoPlasticLinearElastic.h"
#include "include/GooseSolid/PlasticLinearElastic.h"
#include "include/GooseSolid/Miscellaneous.h"

namespace py = pybind11;
namespace GS = GooseSolid;

PYBIND11_PLUGIN(GooseSolid) {

py::module m("GooseSolid","Library with material models");

// -------------------------------------------------------------------------------------------------

py::class_<GS::PlasticLinearElastic>(m,"PlasticLinearElastic")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"   , &GS::PlasticLinearElastic::stress   ,py::arg("eps"))
.def("tangent"  , &GS::PlasticLinearElastic::tangent  ,py::arg("eps"))
.def("increment", &GS::PlasticLinearElastic::increment               )

.def("__repr__",[](const GS::PlasticLinearElastic &a)
  {return "<GooseSolid.PlasticLinearElastic>";});

// -------------------------------------------------------------------------------------------------

py::class_<GS::ViscoPlasticLinearElastic>(m,"ViscoPlasticLinearElastic")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sig0"  ),
  py::arg("gamma0"),
  py::arg("m"     ) = 1.0
)

.def("stress"   , &GS::ViscoPlasticLinearElastic::stress   ,py::arg("eps"),py::arg("dt"))
.def("tangent"  , &GS::ViscoPlasticLinearElastic::tangent  ,py::arg("eps"),py::arg("dt"))
.def("increment", &GS::ViscoPlasticLinearElastic::increment                             )

.def("__repr__",[](const GS::ViscoPlasticLinearElastic &a)
  {return "<GooseSolid.ViscoPlasticLinearElastic>";});

// -------------------------------------------------------------------------------------------------

py::class_<GS::LinearElastic_ViscousFluid>(m,"LinearElastic_ViscousFluid")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy"  ),
  py::arg("Tdamp" ),
  py::arg("Tfluid")
)

.def("stress"     , &GS::LinearElastic_ViscousFluid::stress     ,py::arg("epsdot"),py::arg("dt"))
.def("increment"  , &GS::LinearElastic_ViscousFluid::increment                                  )
.def("setNextSigy", &GS::LinearElastic_ViscousFluid::setNextSigy,py::arg("sigy")                )

.def("__repr__",[](const GS::LinearElastic_ViscousFluid &a)
  {return "<GooseSolid.LinearElastic_ViscousFluid>";});

// -------------------------------------------------------------------------------------------------

m.def("ConvertElasticParameters",&GS::ConvertElasticParameters);
m.def("VonMisesStress"          ,&GS::VonMisesStress          );

// -------------------------------------------------------------------------------------------------

return m.ptr();

} // PYBIND11_PLUGIN
