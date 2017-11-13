
#include <cppmat/cppmat.h>
#include <cppmat/pybind11.h>

#include "../src/GooseMaterial/GooseMaterial.h"

namespace py = pybind11;

PYBIND11_MODULE(GooseMaterial, m)
{

m.doc() = "Library with material models";

// =================================================================================================
// GooseMaterial::Metal
// =================================================================================================

py::module Metal =
  m.def_submodule(
    "Metal", "Material models for metals (etc.)"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain
// =================================================================================================

py::module Metal_LinearStrain =
  Metal.def_submodule(
    "LinearStrain", "Based on the linear strain tensor"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::Elastic
// =================================================================================================

py::module Metal_LinearStrain_Elastic =
  Metal_LinearStrain.def_submodule(
    "Elastic", "Elastic"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::Elastic;

Metal_LinearStrain_Elastic.def("ConvertParameters",&N::ConvertParameters);

}

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d
// =================================================================================================

py::module Metal_LinearStrain_Elastic_Cartesian3d =
  Metal_LinearStrain_Elastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::Elastic::Cartesian3d;

py::class_<N::Material>(Metal_LinearStrain_Elastic_Cartesian3d,"Material")

.def(py::init<double,double>(),
  py::arg("K"     ),
  py::arg("G"     )
)

.def("stress"        , &N::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &N::Material::tangent_stress, py::arg("Eps"))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::NonLinearElastic
// =================================================================================================

py::module Metal_LinearStrain_NonLinearElastic =
  Metal_LinearStrain.def_submodule(
    "NonLinearElastic", "Non-linear elastic"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d
// =================================================================================================

py::module Metal_LinearStrain_NonLinearElastic_Cartesian3d =
  Metal_LinearStrain_NonLinearElastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::NonLinearElastic::Cartesian3d;

py::class_<N::Material>(Metal_LinearStrain_NonLinearElastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double>(),
  py::arg("K"     ),
  py::arg("sig0"  ),
  py::arg("eps0"  ),
  py::arg("n"     )=1.
)

.def("stress"        , &N::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &N::Material::tangent_stress, py::arg("Eps"))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.Elastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoPlastic
// =================================================================================================

py::module Metal_LinearStrain_ElastoPlastic =
  Metal_LinearStrain.def_submodule(
    "ElastoPlastic", "Elasto-plastic"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d
// =================================================================================================

py::module Metal_LinearStrain_ElastoPlastic_Cartesian3d =
  Metal_LinearStrain_ElastoPlastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::ElastoPlastic::Cartesian3d;

py::class_<N::Material>(Metal_LinearStrain_ElastoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &N::Material::stress        , py::arg("Eps"))
.def("tangent_stress", &N::Material::tangent_stress, py::arg("Eps"))
.def("increment"     , &N::Material::increment                     )

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoPlastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic
// =================================================================================================

py::module Metal_LinearStrain_ElastoViscoPlastic =
  Metal_LinearStrain.def_submodule(
    "ElastoViscoPlastic", "Elasto-visco-plastic"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d
// =================================================================================================

py::module Metal_LinearStrain_ElastoViscoPlastic_Cartesian3d =
  Metal_LinearStrain_ElastoViscoPlastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::ElastoViscoPlastic::Cartesian3d;

py::class_<N::Material>(Metal_LinearStrain_ElastoViscoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sig0"  ),
  py::arg("gamma0"),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &N::Material::stress        , py::arg("Eps"), py::arg("dt"))
.def("tangent_stress", &N::Material::tangent_stress, py::arg("Eps"), py::arg("dt"))
.def("increment"     , &N::Material::increment                                    )

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening
// =================================================================================================

py::module Metal_LinearStrain_ElastoViscoPlasticHardening =
  Metal_LinearStrain.def_submodule(
    "ElastoViscoPlasticHardening", "Elasto-visco-plastic, with hardening"
  );

// =================================================================================================
// GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d
// =================================================================================================

py::module Metal_LinearStrain_ElastoViscoPlasticHardening_Cartesian3d =
  Metal_LinearStrain_ElastoViscoPlasticHardening .def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::Metal::LinearStrain::ElastoViscoPlasticHardening::Cartesian3d;

py::class_<N::Material>(Metal_LinearStrain_ElastoViscoPlasticHardening_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("gamma0"),
  py::arg("n"     ),
  py::arg("sigy0" ),
  py::arg("H"     ),
  py::arg("m"     ) = 1.0
)

.def("stress"        , &N::Material::stress        , py::arg("Eps"), py::arg("dt"))
.def("tangent_stress", &N::Material::tangent_stress, py::arg("Eps"), py::arg("dt"))
.def("increment"     , &N::Material::increment                                    )

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.Metal.LinearStrain.ElastoViscoPlasticHardening.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid
// =================================================================================================

py::module AmorphousSolid =
  m.def_submodule(
    "AmorphousSolid", "Material models for amorphous solids"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain
// =================================================================================================

py::module AmorphousSolid_LinearStrain =
  AmorphousSolid.def_submodule(
    "LinearStrain", "Based on the linear strain tensor"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::Elastic
// =================================================================================================

py::module AmorphousSolid_LinearStrain_Elastic =
  AmorphousSolid_LinearStrain.def_submodule(
    "Elastic", "Elastic"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian2d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_Elastic_Cartesian2d =
  AmorphousSolid_LinearStrain_Elastic.def_submodule(
    "Cartesian2d", "2d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian2d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_Elastic_Cartesian2d,"Material")

.def(py::init<double,double>(),
  py::arg("K"),
  py::arg("G")
)

.def("stress",&N::Material::stress, py::arg("Eps"))
.def("energy",&N::Material::energy, py::arg("Eps"))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian2d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_Elastic_Cartesian3d =
  AmorphousSolid_LinearStrain_Elastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::Elastic::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_Elastic_Cartesian3d,"Material")

.def(py::init<double,double>(),
  py::arg("K"),
  py::arg("G")
)

.def("stress",&N::Material::stress, py::arg("Eps"))
.def("energy",&N::Material::energy, py::arg("Eps"))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.Elastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElasticLiquid =
  AmorphousSolid_LinearStrain.def_submodule(
    "ElasticLiquid", "Elastic-Liquid"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElasticLiquid_Cartesian3d =
  AmorphousSolid_LinearStrain_ElasticLiquid.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElasticLiquid::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElasticLiquid_Cartesian3d,"Material")

.def(py::init<double,double,double,double,double>(),
  py::arg("K"     ),
  py::arg("G"     ),
  py::arg("sigy"  ),
  py::arg("Tdamp" ),
  py::arg("Tfluid")
)

.def("stress"     , &N::Material::stress     , py::arg("epsdot"), py::arg("dt"))
.def("increment"  , &N::Material::increment                                    )
.def("setNextSigy", &N::Material::setNextSigy, py::arg("sigy")                 )

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElasticLiquid.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic =
  AmorphousSolid_LinearStrain.def_submodule(
    "ElastoPlastic", "Elasto-plastic"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2d =
  AmorphousSolid_LinearStrain_ElastoPlastic.def_submodule(
    "Cartesian2d", "2d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian2d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian2d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Cartesian2d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3d =
  AmorphousSolid_LinearStrain_ElastoPlastic.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Cartesian3d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::PlanarShear
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear =
  AmorphousSolid_LinearStrain_ElastoPlastic.def_submodule(
    "PlanarShear", "Plasticity along plane"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::PlanarShear::Cartesian2d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear_Cartesian2d =
  AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear.def_submodule(
    "Cartesian2d", "2d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::PlanarShear::Cartesian2d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear_Cartesian2d,"Material")

.def(py::init<double,double,const cppmat::cartesian2d::vector<double> &,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.PlanarShear.Cartesian2d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::PlanarShear::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear_Cartesian3d =
  AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::PlanarShear::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_PlanarShear_Cartesian3d,"Material")

.def(py::init<double,double,const cppmat::cartesian3d::vector<double> &,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.PlanarShear.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth =
  AmorphousSolid_LinearStrain_ElastoPlastic.def_submodule(
    "Smooth", "Elasto-plastic material with a smooth potential"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::Cartesian2d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_Cartesian2d =
  AmorphousSolid_LinearStrain_ElastoPlastic_Smooth.def_submodule(
    "Cartesian2d", "2d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::Cartesian2d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_Cartesian2d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Smooth.Cartesian2d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_Cartesian3d =
  AmorphousSolid_LinearStrain_ElastoPlastic_Smooth.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_Cartesian3d,"Material")

.def(py::init<double,double,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Smooth.Cartesian3d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::PlanarShear
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear =
  AmorphousSolid_LinearStrain_ElastoPlastic_Smooth.def_submodule(
    "PlanarShear", "Plasticity along plane"
  );

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::PlanarShear::Cartesian2d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear_Cartesian2d =
  AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear.def_submodule(
    "Cartesian2d", "2d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::PlanarShear::Cartesian2d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear_Cartesian2d,"Material")

.def(py::init<double,double,const cppmat::cartesian2d::vector<double> &,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Smooth.PlanarShear.Cartesian2d.Material>";});

}

// =================================================================================================
// GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::PlanarShear::Cartesian3d
// =================================================================================================

py::module AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear_Cartesian3d =
  AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear.def_submodule(
    "Cartesian3d", "3d Cartesian coordinate system"
  );

// -------------------------------------------------------------------------------------------------

{

namespace N = GooseMaterial::AmorphousSolid::LinearStrain::ElastoPlastic::Smooth::PlanarShear::Cartesian3d;

py::class_<N::Material>(AmorphousSolid_LinearStrain_ElastoPlastic_Smooth_PlanarShear_Cartesian3d,"Material")

.def(py::init<double,double,const cppmat::cartesian3d::vector<double> &,const std::vector<double> &,bool>(),
  py::arg("K"           ),
  py::arg("G"           ),
  py::arg("n"           ),
  py::arg("eps_y"       )=std::vector<double>(),
  py::arg("init_elastic")=true
)

.def("stress",&N::Material::stress, py::arg("Eps" ))
.def("energy",&N::Material::energy, py::arg("Eps" ))
.def("find"  ,&N::Material::find  , py::arg("epsd"))
.def("eps_y" ,&N::Material::eps_y , py::arg("i"   ))

.def("__repr__",[](const N::Material &a)
  {return "<GooseMaterial.AmorphousSolid.LinearStrain.ElastoPlastic.Smooth.PlanarShear.Cartesian3d.Material>";});

}

// =================================================================================================

} // PYBIND11_PLUGIN
