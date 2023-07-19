# User visible changes in `ConjugateGradient.jl`

# Version 0.3.0

- Status returned by `ConjugateGradient.solve!` is an enumeration which can be compared
  against an integer or another status.

# Version 0.2.0

- Use `NumOptBase` instead of `LazyAlgebra` for basic operations on variables.

- Status returned by `ConjugateGradient.solve!` is a symbolic constant. Method
  `ConjugateGradient.has_converged` is provided to check whether
  `ConjugateGradient.solve!` has converged.

# Version 0.1.0

Proof of concept and first implementation.
