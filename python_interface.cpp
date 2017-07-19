
#include <cppmat/tensor.h>
#include <cppmat/pybind11_tensor.h>

#include "include/GooseSolid/LinearElastic_ViscousFluid.h"

PYBIND11_PLUGIN(GooseSolid) {

py::module m("GooseSolid","Material model: linear elastic / viscous fluid");

// -------------------------------------------------------------------------------------------------

py::class_<GooseSolid::LinearElastic_ViscousFluid>(m,"model")

.def(py::init<double,double,double,double,double,size_t>(),
  pybind11::arg("K"     ),
  pybind11::arg("mu"    ),
  pybind11::arg("sigy"  ),
  pybind11::arg("Tdamp" ),
  pybind11::arg("Tfluid"),
  pybind11::arg("nd"    )=3
)

.def("setNext_sigy", &GooseSolid::LinearElastic_ViscousFluid::setNext_sigy,
  pybind11::arg("sigy")
)

.def("stress", &GooseSolid::LinearElastic_ViscousFluid::stress,
  pybind11::arg("epsdot"),
  pybind11::arg("dt"    )
)

.def("increment", &GooseSolid::LinearElastic_ViscousFluid::increment)

.def("__repr__",[](const GooseSolid::LinearElastic_ViscousFluid &a)
  {return "<GooseSolid.LinearElastic_ViscousFluid>";});

// -------------------------------------------------------------------------------------------------

return m.ptr();

} // PYBIND11_PLUGIN
