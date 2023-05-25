# Linear conjugate gradient algorithm for Julia

[![Build Status](https://github.com/emmt/ConjugateGradient.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/ConjugateGradient.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/ConjugateGradient.jl?svg=true)](https://ci.appveyor.com/project/emmt/ConjugateGradient-jl)
[![Coverage](https://codecov.io/gh/emmt/ConjugateGradient.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/ConjugateGradient.jl)

This [Julia](https://julialang.org/) package provides an implementation of the (linear)
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
