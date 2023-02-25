# Linear conjugate gradient algorithm for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This [Julia][julia-url] package provides an implementation of the (linear)
conjugate gradient algorithm, possibly with a preconditioner.
`ConjugateGradient` exploits
[`LazyAlgebra`](https://emmt.github.io/LazyAlgebra.jl) framework to be as
general as possible without sacrificing performances.

One of the requirements is to avoid allocating resources on every call so that
the method can be used for real-time applications without the risk of being
interrupted by the garbage collector.  To achieve this, the implemented
algorithm has all workspace arrays stored in a context provided by the caller.
This context can be created once and used as many times as wanted (for solving
problems of the same size and type).

The implemented algorithm is very flexible in the type of variables and
operators of the problem and has a number of possible criteria for stopping the
iterations.


## Installation

The easiest way to install `ConjugateGradient.jl` is via Julia registry
[`EmmtRegistry`](https://github.com/emmt/EmmtRegistry):

```julia
using Pkg
pkg"registry add https://github.com/emmt/EmmtRegistry"
pkg"add ConjugateGradient"
```

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/ConjugateGradient.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/ConjugateGradient.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/ConjugateGradient.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/ConjugateGradient.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/ConjugateGradient-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/ConjugateGradient.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/ConjugateGradient.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/ConjugateGradient.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/ConjugateGradient.jl?branch=master

[julia-url]: https://julialang.org/
[julia-pkgs-url]: https://pkg.julialang.org/
