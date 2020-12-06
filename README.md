# Linera conjugate gradient algorithm for Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This [Julia][julia-url] package provides an implementation of the (linear)
conjugate gradient algorithm, possibly with a preconditioner.
`ConjugateGradient` exploits
[`LazyAlgebra`](https://emmt.github.io/LazyAlgebra.jl) framework to be as
general as possible without sacrificing performances.

## Installation

`ConjugateGradient` is not yet an [official Julia package][julia-pkgs-url] but
it is easy to install.  In Julia, hit the `]` key to switch to the package
manager REPL (you should get a `... pkg>` prompt) and type:

```julia
... pkg> add https://emmt.github.io/LazyAlgebra.jl
... pkg> add https://emmt.github.io/ConjugateGradient.jl
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
