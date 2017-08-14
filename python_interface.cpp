
#include <cppmat/tensor.h>
#include <cppmat/pybind11_tensor.h>

#include "include/GooseSolid/LinearElastic.h"
#include "include/GooseSolid/NonLinearElastic.h"
#include "include/GooseSolid/PlasticLinearElastic.h"
#include "include/GooseSolid/ViscoPlasticLinearElastic.h"
#include "include/GooseSolid/ViscoPlasticHardeningLinearElastic.h"
#include "include/GooseSolid/LinearElastic_ViscousFluid.h"
#include "include/GooseSolid/Miscellaneous.h"

namespace py = pybind11;
namespace GS = GooseSolid;

PYBIND11_PLUGIN(GooseSolid) {

py::module m("GooseSolid","Library with material models");

// -------------------------------------------------------------------------------------------------

py::class_<GS::LinearElastic>(m,"LinearElastic")

.def(py::init<double,double>(),
  py::arg("K"     ),
  py::arg("G"     )
)

.def("stress"        , &GS::LinearElastic::stress        , py::arg("eps"))
.def("tangent_stress", &GS::LinearElastic::tangent_stress, py::arg("eps"))
.def("tangent"       , &GS::LinearElastic::tangent       , py::arg("eps"))

.def("__repr__",[](const GS::LinearElastic &a)
  {return "<GooseSolid.LinearElastic>";});

// -------------------------------------------------------------------------------------------------

py::class_<GS::NonLinearElastic>(m,"NonLinearElastic")

.def(py::init<double,double,double,double>(),
  py::arg("K"     ),
  py::arg("sig0"  ),
  py::arg("eps0"  ),
  py::arg("n"     )=1.
)

.def("stress"        , &GS::NonLinearElastic::stress        , py::arg("eps"))
.def("tangent_stress", &GS::NonLinearElastic::tangent_stress, py::arg("eps"))
.def("tangent"       , &GS::NonLinearElastic::tangent       , py::arg("eps"))

.def("__repr__",[](const GS::NonLinearElastic &a)
  {return "<GooseSolid.NonLinearElastic>";});

// -------------------------------------------------------------------------------------------------

py::class_<GS::PlasticLinearElastic>(m,"PlasticLinearElastic")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GS::PlasticLinearElastic::stress        , py::arg("eps"))
.def("tangent_stress", &GS::PlasticLinearElastic::tangent_stress, py::arg("eps"))
.def("tangent"       , &GS::PlasticLinearElastic::tangent       , py::arg("eps"))
.def("increment"     , &GS::PlasticLinearElastic::increment                     )

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

.def("stress"        , &GS::ViscoPlasticLinearElastic::stress        , py::arg("eps"), py::arg("dt"))
.def("tangent_stress", &GS::ViscoPlasticLinearElastic::tangent_stress, py::arg("eps"), py::arg("dt"))
.def("tangent"       , &GS::ViscoPlasticLinearElastic::tangent       , py::arg("eps"), py::arg("dt"))
.def("increment"     , &GS::ViscoPlasticLinearElastic::increment                                    )

.def("__repr__",[](const GS::ViscoPlasticLinearElastic &a)
  {return "<GooseSolid.ViscoPlasticLinearElastic>";});

// -------------------------------------------------------------------------------------------------

py::class_<GS::ViscoPlasticHardeningLinearElastic>(m,"ViscoPlasticHardeningLinearElastic")

.def(py::init<double,double,double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("gamma0"),
  py::arg("n"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GS::ViscoPlasticHardeningLinearElastic::stress        , py::arg("eps"), py::arg("dt"))
.def("tangent_stress", &GS::ViscoPlasticHardeningLinearElastic::tangent_stress, py::arg("eps"), py::arg("dt"))
.def("tangent"       , &GS::ViscoPlasticHardeningLinearElastic::tangent       , py::arg("eps"), py::arg("dt"))
.def("increment"     , &GS::ViscoPlasticHardeningLinearElastic::increment                                    )

.def("__repr__",[](const GS::ViscoPlasticHardeningLinearElastic &a)
  {return "<GooseSolid.ViscoPlasticHardeningLinearElastic>";});


// -------------------------------------------------------------------------------------------------

py::class_<GS::LinearElastic_ViscousFluid>(m,"LinearElastic_ViscousFluid")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy"  ),
  py::arg("Tdamp" ),
  py::arg("Tfluid")
)

.def("stress"     , &GS::LinearElastic_ViscousFluid::stress     , py::arg("epsdot"), py::arg("dt"))
.def("increment"  , &GS::LinearElastic_ViscousFluid::increment                                    )
.def("setNextSigy", &GS::LinearElastic_ViscousFluid::setNextSigy, py::arg("sigy")                 )

.def("__repr__",[](const GS::LinearElastic_ViscousFluid &a)
  {return "<GooseSolid.LinearElastic_ViscousFluid>";});

// -------------------------------------------------------------------------------------------------

m.def("ConvertElasticParameters",&GS::ConvertElasticParameters);

// -------------------------------------------------------------------------------------------------

return m.ptr();

} // PYBIND11_PLUGIN
