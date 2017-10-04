
#include <cppmat/tensor2.h>
#include <cppmat/tensor3.h>
#include <cppmat/pybind11_tensor2.h>
#include <cppmat/pybind11_tensor3.h>

#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElasticLiquid/Cartesian3d.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian3d.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian3dPlanarShear.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2dPlanarShear.h"
#include "../src/GooseMaterial/Metal/LinearStrain/Elastic/Cartesian3d.h"
#include "../src/GooseMaterial/Metal/LinearStrain/Elastic/miscellaneous.h"
#include "../src/GooseMaterial/Metal/LinearStrain/ElastoPlastic/Cartesian3d.h"
#include "../src/GooseMaterial/Metal/LinearStrain/ElastoViscoPlastic/Cartesian3d.h"
#include "../src/GooseMaterial/Metal/LinearStrain/ElastoViscoPlasticHardening/Cartesian3d.h"
#include "../src/GooseMaterial/Metal/LinearStrain/NonLinearElastic/Cartesian3d.h"

namespace py = pybind11;

PYBIND11_MODULE(GooseMaterial, m)
{

m.doc() = "Library with material models";

// GooseMaterial::*
py::module Metal                                                            = m                                              .def_submodule( "Metal"                       , "Material models for metals (etc.)"                      );
py::module AmorphousSolid                                                   = m                                              .def_submodule( "AmorphousSolid"              , "Material models for amorphous solids"                   );
// GooseMaterial::*::*
py::module Metal_LinearStrain                                               = Metal                                          .def_submodule( "LinearStrain"                , "Based on the linear strain tensor"                      );
py::module AmorphousSolid_LinearStrain                                      = AmorphousSolid                                 .def_submodule( "LinearStrain"                , "Based on the linear strain tensor"                      );
// GooseMaterial::*::*::*
py::module Metal_LinearStrain_Elastic                                       = Metal_LinearStrain                             .def_submodule( "Elastic"                     , "Elastic"                                                );
py::module Metal_LinearStrain_NonLinearElastic                              = Metal_LinearStrain                             .def_submodule( "NonLinearElastic"            , "Non-linear elastic"                                     );
py::module Metal_LinearStrain_ElastoPlastic                                 = Metal_LinearStrain                             .def_submodule( "ElastoPlastic"               , "Elasto-plastic"                                         );
py::module Metal_LinearStrain_ElastoViscoPlastic                            = Metal_LinearStrain                             .def_submodule( "ElastoViscoPlastic"          , "Elasto-visco-plastic"                                   );
py::module Metal_LinearStrain_ElastoViscoPlasticHardening                   = Metal_LinearStrain                             .def_submodule( "ElastoViscoPlasticHardening" , "Elasto-visco-plastic, with hardening"                   );
py::module AmorphousSolid_LinearStrain_ElasticLiquid                        = AmorphousSolid_LinearStrain                    .def_submodule( "ElasticLiquid"               , "Elastic-Liquid"                                         );
py::module AmorphousSolid_LinearStrain_ElastoPlastic                        = AmorphousSolid_LinearStrain                    .def_submodule( "ElastoPlastic"               , "Elasto-plastic"                                         );
// GooseMaterial::*::*::*::*
py::module Metal_LinearStrain_Elastic_Cartesian3d                           = Metal_LinearStrain_Elastic                     .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module Metal_LinearStrain_NonLinearElastic_Cartesian3d                  = Metal_LinearStrain_NonLinearElastic            .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module Metal_LinearStrain_ElastoPlastic_Cartesian3d                     = Metal_LinearStrain_ElastoPlastic               .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module Metal_LinearStrain_ElastoViscoPlastic_Cartesian3d                = Metal_LinearStrain_ElastoViscoPlastic          .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module Metal_LinearStrain_ElastoViscoPlasticHardening_Cartesian3d       = Metal_LinearStrain_ElastoViscoPlasticHardening .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module AmorphousSolid_LinearStrain_ElasticLiquid_Cartesian3d            = AmorphousSolid_LinearStrain_ElasticLiquid      .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3d            = AmorphousSolid_LinearStrain_ElastoPlastic      .def_submodule( "Cartesian3d"                 , "3d Cartesian coordinate system"                         );
py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2d            = AmorphousSolid_LinearStrain_ElastoPlastic      .def_submodule( "Cartesian2d"                 , "2d Cartesian coordinate system"                         );
py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3dPlanarShear = AmorphousSolid_LinearStrain_ElastoPlastic      .def_submodule( "Cartesian3dPlanarShear"      , "3d Cartesian coordinate system - plasticity along plane");
py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2dPlanarShear = AmorphousSolid_LinearStrain_ElastoPlastic      .def_submodule( "Cartesian2dPlanarShear"      , "2d Cartesian coordinate system - plasticity along plane");

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material>(Metal_LinearStrain_Elastic_Cartesian3d,"Material")

.def(py::init<double,double>(),
  py::arg("K"     ),
  py::arg("G"     )
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material::tangent_stress, py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

// -------------------------------------------------------------------------------------------------

Metal_LinearStrain_Elastic.def("ConvertParameters",&GooseMaterial::Metal::LinearStrain::Elastic::ConvertParameters);

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material>(Metal_LinearStrain_NonLinearElastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double>(),
  py::arg("K"     ),
  py::arg("sig0"  ),
  py::arg("eps0"  ),
  py::arg("n"     )=1.
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material::tangent_stress, py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material>(Metal_LinearStrain_ElastoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::tangent_stress, py::arg("Eps"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::increment                     )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoPlastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material>(Metal_LinearStrain_ElastoViscoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sig0"  ),
  py::arg("gamma0"),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::stress        , py::arg("Eps"), py::arg("dt"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::tangent_stress, py::arg("Eps"), py::arg("dt"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::increment                                    )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material>(Metal_LinearStrain_ElastoViscoPlasticHardening_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("gamma0"),
  py::arg("n"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::stress        , py::arg("Eps"), py::arg("dt"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::tangent_stress, py::arg("Eps"), py::arg("dt"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::increment                                    )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlasticHardening.Cartesian3d.Material>";});


// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material>(AmorphousSolid_LinearStrain_ElasticLiquid_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy"  ),
  py::arg("Tdamp" ),
  py::arg("Tfluid")
)

.def("stress"     , &GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material::stress     , py::arg("epsdot"), py::arg("dt"))
.def("increment"  , &GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material::increment                                    )
.def("setNextSigy", &GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material::setNextSigy, py::arg("sigy")                 )

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElasticLiquid.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress"  ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::stress   , py::arg("Eps"))
.def("eps_d"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::eps_d    , py::arg("Eps"))
.def("eps_m"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::eps_m    , py::arg("Eps"))
.def("energy"  , py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy  ), py::arg("Eps"))
.def("energy_m", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy_m), py::arg("Eps"))
.def("energy_d", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy_d), py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress"  ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::stress   , py::arg("Eps"))
.def("eps_d"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::eps_d    , py::arg("Eps"))
.def("eps_m"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::eps_m    , py::arg("Eps"))
.def("energy"  , py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy  ), py::arg("Eps"))
.def("energy_m", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy_m), py::arg("Eps"))
.def("energy_d", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy_d), py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian2d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3dPlanarShear,"Material")

.def(py::init<double,double,const cppmat::vector3<double>&,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress"  ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::stress   , py::arg("Eps"))
.def("eps_d"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::eps_d    , py::arg("Eps"))
.def("eps_s"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::eps_s    , py::arg("Eps"))
.def("eps_n"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::eps_n    , py::arg("Eps"))
.def("eps_m"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::eps_m    , py::arg("Eps"))
.def("energy"  , py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::energy  ), py::arg("Eps"))
.def("energy_m", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::energy_m), py::arg("Eps"))
.def("energy_s", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::energy_s), py::arg("Eps"))
.def("energy_n", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material::energy_n), py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3dPlanarShear::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian3dPlanarShear.Material>";});

// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2dPlanarShear,"Material")

.def(py::init<double,double,const cppmat::vector2<double>&,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress"  ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::stress   , py::arg("Eps"))
.def("eps_d"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::eps_d    , py::arg("Eps"))
.def("eps_s"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::eps_s    , py::arg("Eps"))
.def("eps_n"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::eps_n    , py::arg("Eps"))
.def("eps_m"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::eps_m    , py::arg("Eps"))
.def("energy"  , py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::energy  ), py::arg("Eps"))
.def("energy_m", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::energy_m), py::arg("Eps"))
.def("energy_s", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::energy_s), py::arg("Eps"))
.def("energy_n", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material::energy_n), py::arg("Eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2dPlanarShear::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian2dPlanarShear.Material>";});

// =================================================================================================

} // PYBIND11_PLUGIN
