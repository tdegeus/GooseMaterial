
#include <cppmat/tensor2.h>
#include <cppmat/tensor3.h>
#include <cppmat/pybind11_tensor2.h>
#include <cppmat/pybind11_tensor3.h>

#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElasticLiquid/Cartesian3d.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian3d.h"
#include "../src/GooseMaterial/AmorphousSolid/LinearStrain/ElastoPlastic/Cartesian2d.h"
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

py::module Metal = m.def_submodule(
  "Metal",
  "Material models for metal (etc.)"
);

py::module MetalLinearStrain = Metal.def_submodule(
  "LinearStrain",
  "Material models based on the linear strain tensor"
);

py::module MetalLinearStrainElastic = MetalLinearStrain.def_submodule(
  "Elastic",
  "Elastic material model"
);

py::module MetalLinearStrainElasticCartesian3d = MetalLinearStrainElastic.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module MetalLinearStrainNonLinearElastic = MetalLinearStrain.def_submodule(
  "NonLinearElastic",
  "Non-linear elastic material model"
);

py::module MetalLinearStrainNonLinearElasticCartesian3d = MetalLinearStrainNonLinearElastic.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module MetalLinearStrainElastoPlastic = MetalLinearStrain.def_submodule(
  "ElastoPlastic",
  "Elasto-plastic material model"
);

py::module MetalLinearStrainElastoPlasticCartesian3d = MetalLinearStrainElastoPlastic.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module MetalLinearStrainElastoViscoPlastic = MetalLinearStrain.def_submodule(
  "ElastoViscoPlastic",
  "Elasto-visco-plastic material model"
);

py::module MetalLinearStrainElastoViscoPlasticCartesian3d = MetalLinearStrainElastoViscoPlastic.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module MetalLinearStrainElastoViscoPlasticHardening = MetalLinearStrain.def_submodule(
  "ElastoViscoPlasticHardening",
  "Elasto-visco-plastic material model, with hardening"
);

py::module MetalLinearStrainElastoViscoPlasticHardeningCartesian3d = MetalLinearStrainElastoViscoPlasticHardening.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module AmorphousSolid = m.def_submodule(
  "AmorphousSolid",
  "Material models for amorphous solids"
);

py::module AmorphousSolidLinearStrain = AmorphousSolid.def_submodule(
  "LinearStrain",
  "Material models based on the linear strain tensor"
);

py::module AmorphousSolidLinearStrainElasticLiquid = AmorphousSolidLinearStrain.def_submodule(
  "ElasticLiquid",
  "Elastic-Liquid material model"
);

py::module AmorphousSolidLinearStrainElasticLiquidCartesian3d = AmorphousSolidLinearStrainElasticLiquid.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module AmorphousSolidLinearStrainElastoPlastic = AmorphousSolidLinearStrain.def_submodule(
  "ElastoPlastic",
  "Elasto-plastic material model"
);

py::module AmorphousSolidLinearStrainElastoPlasticCartesian3d = AmorphousSolidLinearStrainElastoPlastic.def_submodule(
  "Cartesian3d",
  "Defined on a 3d Cartesian coordinate system"
);

py::module AmorphousSolidLinearStrainElastoPlasticCartesian2d = AmorphousSolidLinearStrainElastoPlastic.def_submodule(
  "Cartesian2d",
  "Defined on a 3d Cartesian coordinate system"
);

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material>(MetalLinearStrainElasticCartesian3d,"Material")

.def(py::init<double,double>(),
  py::arg("K"     ),
  py::arg("G"     )
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material::stress        , py::arg("eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material::tangent_stress, py::arg("eps"))

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

// -------------------------------------------------------------------------------------------------

MetalLinearStrainElastic.def("ConvertParameters",&GooseMaterial::Metal::LinearStrain::Elastic::ConvertParameters);

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material>(MetalLinearStrainNonLinearElasticCartesian3d,"Material")

.def(py::init<double,double,double,double>(),
  py::arg("K"     ),
  py::arg("sig0"  ),
  py::arg("eps0"  ),
  py::arg("n"     )=1.
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material::stress        , py::arg("eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material::tangent_stress, py::arg("eps"))

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material>(MetalLinearStrainElastoPlasticCartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::stress        , py::arg("eps"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::tangent_stress, py::arg("eps"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material::increment                     )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoPlastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material>(MetalLinearStrainElastoViscoPlasticCartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sig0"  ),
  py::arg("gamma0"),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::stress        , py::arg("eps"), py::arg("dt"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::tangent_stress, py::arg("eps"), py::arg("dt"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material::increment                                    )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material>(MetalLinearStrainElastoViscoPlasticHardeningCartesian3d,"Material")

.def(py::init<double,double,double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("gamma0"),
  py::arg("n"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::stress        , py::arg("eps"), py::arg("dt"))
.def("tangent_stress", &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::tangent_stress, py::arg("eps"), py::arg("dt"))
.def("increment"     , &GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material::increment                                    )

.def("__repr__",[](const GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlasticHardening.Cartesian3d.Material>";});


// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d::Material>(AmorphousSolidLinearStrainElasticLiquidCartesian3d,"Material")

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

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material>(AmorphousSolidLinearStrainElastoPlasticCartesian3d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       ),
  py::arg("init_elastic")=true
)

.def("stress"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::stress    , py::arg("eps"))
.def("eps_eq"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::eps_eq    , py::arg("eps"))
.def("eps_m"    ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::eps_m     , py::arg("eps"))
.def("energy"   , py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy   ), py::arg("eps"))
.def("energy_m ", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy_m ), py::arg("eps"))
.def("energy_eq", py::overload_cast<const cppmat::tensor3_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material::energy_eq), py::arg("eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian3d.Material>";});

// =================================================================================================

py::class_<GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material>(AmorphousSolidLinearStrainElastoPlasticCartesian2d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       ),
  py::arg("init_elastic")=true
)

.def("stress"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::stress    , py::arg("eps"))
.def("eps_eq"   ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::eps_eq    , py::arg("eps"))
.def("eps_m"    ,                                                      &GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::eps_m     , py::arg("eps"))
.def("energy"   , py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy   ), py::arg("eps"))
.def("energy_m ", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy_m ), py::arg("eps"))
.def("energy_eq", py::overload_cast<const cppmat::tensor2_2s<double>&>(&GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material::energy_eq), py::arg("eps"))

.def("__repr__",[](const GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian2d.Material>";});

// =================================================================================================

} // PYBIND11_PLUGIN
