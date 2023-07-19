module ConjugateGradientTests

using Test, LinearAlgebra, ConjugateGradient
using NumOptBase: Diag, norm2

VERBOSE = true
io = (VERBOSE ? stdout : ConjugateGradient.quiet)
dims = (29,11)
nrows, ncols = dims
status_successful = (ConjugateGradient.F_TEST_SATISFIED,
                     ConjugateGradient.G_TEST_SATISFIED,
                     ConjugateGradient.X_TEST_SATISFIED)
status_ok = (status_successful...,
             ConjugateGradient.TOO_MANY_ITERATIONS)

@testset "Status" begin
    for x in instances(ConjugateGradient.Status)
        @test ConjugateGradient.has_converged(x) == (x ∈ status_successful)
        y = @inferred Symbol(x)
        @test ConjugateGradient.Status(y) === x
    end
end

@testset "least-square fit ($T)" for T in (Float32, Float64)
    # State the problem (least square fit).
    H = rand(T, dims) # direct model
    x = rand(T, ncols) # ground truth
    y = H*x + T(0.01)*rand(T, nrows) # noisy data

    # LHS matrix and RHS vector of the normal equations.
    A = H'*H
    b = H'*y

    # Least-squares solution using Julia linear algebra.
    x0 = A\b

    # Tolerances.
    rtol = 1e-3
    atol = rtol*norm2(x0)

    # Solve with conjugate gradient algorithm.
    n = length(b)
    ws = ConjugateGradient.Context(b;
                                   precond = true,
                                   maxiter = n+5,
                                   ftol = (0.0, 0.0),
                                   gtol = (0.0, 0.0),
                                   xtol = (0.0, 0.0))
    x1 = fill!(similar(b), 0)
    VERBOSE && println("\n# Run the conjugate gradient from 0:")
    status = ConjugateGradient.solve!(x1, A, b, ws, io)
    @test status ∈ status_ok
    VERBOSE && println("rel. err. for x1: ", norm2(x1 - x0)/norm2(x0))
    @test norm2(x1 - x0) ≤ atol

    # Simple diagonal preconditioner.
    q = convert(Vector{T}, diag(A))
    for i in eachindex(q)
        q[i] = 1/q[i]
    end
    M = Diag(q)

    x2 = fill!(similar(b), 0)
    VERBOSE && println("\n# Run the preconditioned conjugate gradient from 0:")
    status = ConjugateGradient.solve!(x2, A, b, M, ws, io)
    @test status ∈ status_ok
    VERBOSE && println("rel. err. for x1: ", norm2(x1 - x0)/norm2(x0))
    @test norm2(x1 - x0) ≤ atol
end

end # module
