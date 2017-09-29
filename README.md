# GooseSolid

This library includes several material models of solids, to be used for example in FEM programs. The idea is that each material model is included in a separate header-only C++ library, which is completely independent. As a rule of principle the different material model follow the same structure and use [<cppmat/tensor.h>](https://github.com/tdegeus/cppmat) and/or [<cppmat/tensor3.h>](https://github.com/tdegeus/cppmat) as tensor library, but the details may vary from material to material. 

**Contents**

<!-- MarkdownTOC -->

- [Materials](#materials)
- [Installation](#installation)
- [C++](#c)
- [Python](#python)

<!-- /MarkdownTOC -->

>   **Disclaimer**
>   
>   This library is free to use under the [GPLv3 license](https://github.com/tdegeus/GooseSolid/blob/master/LICENSE). Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bug reports or feature requests can be filed on [GitHub](https://github.com/tdegeus/GooseSolid). As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.
>   
>   Download: [.zip file](https://github.com/tdegeus/GooseSolid/zipball/master) | [.tar.gz file](https://github.com/tdegeus/GooseSolid/tarball/master).
>   
>   (c - [GPLv3](https://github.com/tdegeus/GooseSolid/blob/master/LICENSE)) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | [github.com/tdegeus/GooseSolid](https://github.com/tdegeus/GooseSolid)

## Materials

The following materials have been implemented:

*   [GooseSolid::LinearElastic](./docs/LinearElastic/LinearElastic.pdf)
*   [GooseSolid::NonLinearElastic](./docs/NonLinearElastic/NonLinearElastic.pdf)
*   [GooseSolid::PlasticLinearElastic](./docs/PlasticLinearElastic/PlasticLinearElastic.pdf)
*   [GooseSolid::ViscoPlasticHardeningLinearElastic](./docs/ViscoPlasticHardeningLinearElastic/ViscoPlasticHardeningLinearElastic.pdf)
*   [GooseSolid::ViscoPlasticLinearElastic](./docs/ViscoPlasticLinearElastic/ViscoPlasticLinearElastic.pdf)
*   [GooseSolid::LinearElastic_ViscousFluid](./docs/LinearElastic_ViscousFluid/LinearElastic_ViscousFluid.pdf)

## Installation

This library is header only, so in principle one does not need to install anything to be able to use it in C++. However it might be a good idea to copy the header to the appropriate location on the system, and to configure `pkg-config`. The steps to follow are identical to those of [cppmat](http://cppmat.geus.me/en/latest/compile.html). The most common thing to do is:

```bash
cd /path/to/GooseSolid
mkdir src
mkdir build
cd build
cmake ..
make install
```

## C++

The use of the header-only C++ library can be illustrated using an example:

```cppmat
#include <GooseSolid/LinearElastic.h>

include main()
{
  ...

  // elastic parameters
  double E  = 1.0;
  double nu = 0.3;
  double K,G;

  // convert elastic parameters to the pair required by "GooseSolid::LinearElastic"
  std::tie(K,G) = GooseSolid::ConvertElasticParameters("E,nu",E,nu,"K,G");

  // linear elastic material
  GooseSolid::LinearElastic mat = GooseSolid::LinearElastic(K,G);

  ...
}
```

Whereby in the case of an inhomogeneous solid one wants to have many material definitions. For example one per finite element:

```cpp
#include <vector>
#include <GooseSolid/LinearElastic.h>

int main()
{

  ...

  std::vector<GooseSolid::LinearElastic> mat( nelem );

  for ( size_t i = 0 ; i < nelem ; ++i )
  {
    double K = ...
    double G = ...

    GooseSolid::LinearElastic emat = GooseSolid::LinearElastic(K,G);

    mat.push_back( emat );
  }

  ...
}
```

Compilation is easy, since the library is header only. Provided that ``pkg-config`` is setup on can use:

```bash
clang++ `pkg-config --cflags cppmat GooseSolid` -std=c++14 -o example example_name.cpp
```

Or one can use

```bash
clang++ -I/path/to/cppmat -I/path/to/GooseSolid -std=c++14 -o example example_name.cpp
```

## Python

The materials also have a Python interface. To use it just run

```bash
python3 setup.py build
python3 setup.py install
```

where one has to replace the executable ``python3`` with ones favorite Python version (e.g. ``python``).


