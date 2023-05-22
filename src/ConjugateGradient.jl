module ConjugateGradient

using LazyAlgebra
using LinearAlgebra
using Printf

"""
    ConjugateGradient.reason(status) -> str

yields a textual description of the status returned by the conjugate gradient
method.

"""
reason(status::Symbol) =
    status === :NOT_POSITIVE_DEFINITE ? "LHS operator is not positive definite" :
    status === :TOO_MANY_ITERATIONS   ? "too many iterations" :
    status === :F_TEST_SATISFIED      ? "function reduction test satisfied" :
    status === :G_TEST_SATISFIED      ? "gradient test satisfied" :
    status === :X_TEST_SATISFIED      ? "variables change test satisfied" :
    "unknown conjugate gradient result"

"""
    ConjugateGradient.has_converged(status) -> bool

yields whether ConjugateGradient.solve!` has converged according to `status`.
method.

"""
has_converged(status::Symbol) =
    (status === :F_TEST_SATISFIED) | (status === :G_TEST_SATISFIED) | (status === :X_TEST_SATISFIED)

struct Context{V}
    p::V # current search direction
    q::V # q = A*p
    r::V # current residuals: r = A*x - b
    z::V # preconditioned residuals: z = M*r
    maxiter::Int
    restart::Int
    fatol::Float64 # absolute tolerance for function reduction
    frtol::Float64 # relative tolerance for function reduction
    gatol::Float64 # absolute tolerance for gradient norm
    grtol::Float64 # relative tolerance for gradient norm
    xatol::Float64 # absolute tolerance for norm of variables change
    xrtol::Float64 # absolute tolerance for norm of variables change
    function Context{V}(p::V, q::V, r::V, z::V=r;
                        maxiter::Integer = typemax(Int),
                        restart::Integer = min(50, length(p)),
                        ftol::NTuple{2,Real} = (0.0, 1e-8),
                        gtol::NTuple{2,Real} = (0.0, 1e-5),
                        xtol::NTuple{2,Real} = (0.0, 1e-6)) where {V}
        maxiter ≥ 0 || bad_argument(
            "bad maximum number of iterations (maxiter = ", maxiter, ")")
        restart ≥ 1 || bad_argument(
            "bad number of iterations for restarting (restart = ", restart,")")
        ftol[1] ≥ 0 || bad_argument(
            "bad function reduction absolute tolerance (ftol[1] = ", ftol[1], ")")
        0 ≤ ftol[2] < 1 || bad_argument(
            "bad function reduction relative tolerance (ftol[2] = ", ftol[2], ")")
        gtol[1] ≥ 0 || bad_argument(
            "bad gradient absolute tolerance (gtol[1] = ", gtol[1], ")")
        0 ≤ gtol[2] < 1 || bad_argument(
            "bad gradient relative tolerance (gtol[2] = ", gtol[2], ")")
        xtol[1] ≥ 0 || bad_argument(
            "bad variables change absolute tolerance (xtol[1] = ", xtol[1], ")")
        0 ≤ xtol[2] < 1 || bad_argument(
            "bad variables change relative tolerance (xtol[2] = ", xtol[1], ")")
        return new{V}(p, q, r, z, maxiter, restart, ftol[1], ftol[2],
                      gtol[1], gtol[2], xtol[1], xtol[2])
    end
end

function Base.show(io::IO, ctx::Context{V}) where {V}
    print(io, "ConjugateGradient.Context(")
    print(io, "x::$V; # array of size $(size(ctx.r))\n")
    print(io, "    precond = $(ctx.z !== ctx.r),\n")
    print(io, "    maxiter = $(ctx.maxiter),\n")
    print(io, "    restart = $(ctx.restart),\n")
    print(io, "    ftol = ($(ctx.fatol), $(ctx.frtol)),\n")
    print(io, "    gtol = ($(ctx.gatol), $(ctx.grtol)),\n")
    print(io, "    xtol = ($(ctx.xatol), $(ctx.xrtol)))")
end

"""
    ConjugateGradient.Context(x) -> ctx

yields a structure with all parameters and storage for temporary variables
needed for running the linear conjugate gradient algorithm. Argument `x`
specifies the variables of the problem and is used to allocate temporary
variables by calling `LazyAlgebra.vcreate(x)`. The returned context holds no
references on `x`. Algorithm parameters are specified by the following
keywords:

* `precond` is to specify whether to allocate temporary variables to store the
  preconditioned residuals. By default, `precond = false`. If false, only the
  un-preconditioned version of the algorithm can be run.

* `maxiter` is to specify the maximum number of iterations to perform which is
  practically unlimited by default.

* `restart` is to specify the number of consecutive iterations before
  restarting the conjugate gradient recurrence. Restarting the algorithm is to
  cope with the accumulation of rounding errors. By default, `restart =
  min(50,length(x)+1)`. Set `restart` to a value less or equal zero or greater
  than `maxiter` if you do not want that any restarts ever occur.

* `ftol = (fatol,frtol)` is to specify the absolute and relative tolerances for
  the function reduction. By default, `ftol = (0.0,1e-8)`.

* `gtol = (gatol,grtol)` is to specify the absolute and relative tolerances for
  stopping the algorithm based on the gradient of the objective function.
  Convergence occurs when the Mahalanobis norm of the residuals (which is that
  of the gradient of the associated objective function) is less or equal the
  largest of `gatol` and `grtol` times the Mahalanobis norm of the initial
  residuals. By default, `gtol = (0.0,1e-5)`.

* `xtol = (xatol,xrtol)` is to specify the absolute and relative tolerances for
  the change in variables. By default, `xtol = (0.0,1e-6)`.

Just recall the constructor with an instance of `ConjugateGradient.Context` to
re-use the same temporary variables (if possible) but different parameters:

    ConjugateGradient.Context(ctx; precond=…, maxiter=…, ...)

"""
function Context(x::V; precond::Bool = false, kwds...) where {V}
    p = vcreate(x)
    q = vcreate(x)
    r = vcreate(x)
    z = (precond ? vcreate(x) : r)
    return Context{V}(p, q, r, z; kwds...)
end

function Context(ctx::Context{V};
                 precond::Bool = false,
                 maxiter::Integer = ctx.maxiter,
                 restart::Integer = ctx.restart,
                 ftol::NTuple{2,Real} = (ctx.fatol, ctx.frtol),
                 gtol::NTuple{2,Real} = (ctx.gatol, ctx.grtol),
                 xtol::NTuple{2,Real} = (ctx.xatol, ctx.xrtol)) where {V}
    return Context{V}(ctx.p,
                      ctx.q,
                      ctx.r,
                      (precond && ctx.z === ctx.r ? vcreate(ctx.r) : ctx.z);
                      maxiter = maxiter,
                      restart = restart,
                      ftol = ftol,
                      gtol = gtol,
                      xtol = xtol)
end

"""
    ConjugateGradient.solve!(x, A, b[, M], ctx[, io=devnull]) -> status

runs the (preconditioned) linear conjugate gradient algorithm to solve the
system of equations `A*x = b` in `x` and according to the settings in `ctx`.

Argument `x` stores the initial solution on entry and the estimated solution on
return.

Argument `A` implements the *left-hand-side (LHS) matrix* of the equations. It
is used as `LazyAlgebra.vmul!(dst,A,src)` to store in `dst` the result of
applying `A` to `src` and where `src` and `dst` are similar to arguments `x`
and `b`. If none of these is suitable, the method `LazyAlgebra.vmul!` can be
extended. Note that, as `A` and `M` must be symmetric, it may be faster to
apply their adjoint.

Argument `b` is the *right-hand-side (RHS) vector* of the equations.  It is
left unchanged.

Argument `ctx` is a [`ConjugateGradient.Context`](@ref) structure storing all
temporary variables and parameters of the algorithm. This argument is reusable
and is required to avoid any additional allocations.

Optional argument `M` is a preconditioner. If `M` is unspecified or if `M` is
`LazyAlgebra.Identity()`, the unpreconditioned version of the algorithm is run.
The preconditioner can be specified in various forms (as for the LHS operator
`A`).

Optional argument `io` is to specify an `IO` instance to which print various
information at each iterations. Nothing is printed if `io` is `devnull` which
is the default.


## Convergence criteria

Provided `A` be positive definite, the solution `x` of the equations `A*x = b`
is unique and is also the minimum of the following convex quadratic objective
function:

    f(x) = (1/2)*x'*A*x - b'*x + ϵ

where `ϵ` is an arbitrary constant. The gradient of this objective function is:

    ∇f(x) = A*x - b

hence solving `A*x = b` for `x` yields the minimum of `f(x)`. The variations of
`f(x)` between successive iterations, the norm of the gradient `∇f(x)`, or the
norm of the variation of variables `x` may be used to decide the convergence of
the algorithm (see keywords `ftol`, `gtol` and `xtol` below).

Let `x_{k}`, `f_{k} = f(x_{k})` and `∇f_{k} = ∇f(x_{k})` denote the variables,
the objective function and its gradient at iteration `k`. The argument `x`
gives the initial variables `x_{0}`. Starting with `k = 0`, the different
possibilities for the convergence of the algorithm are listed below.

* The convergence in the function reduction between succesive iterations occurs
  at iteration `k ≥ 1` if:

  ```
  f_{k-1} - f_{k} ≤ max(fatol, frtol*max_{k' ≤ k}(f_{k'-1} - f_{k'}))
  ```

* The convergence in the gradient norm occurs at iteration `k ≥ 0` if:

  ```
  ‖∇f_{k}‖_M ≤ max(gatol, grtol*‖∇f_{0}‖_M)
  ```

  where `‖u‖_M = sqrt(u'*M*u)` is the Mahalanobis norm of `u` with precision
  matrix `M` which is equal to the usual Euclidean norm of `u` if no
  preconditioner is used or if `M` is the identity.

* The convergence in the variables occurs at iteration `k ≥ 1` if:

  ```
  ‖x_{k} - x_{k-1}‖ ≤ max(xatol, xrtol*‖x_{k}‖)
  ```

In the conjugate gradient algorithm, the objective function is always reduced
at each iteration, but be aware that the gradient and the change of variables
norms are not always reduced at each iteration.


## Returned Status

The returned value `status` is one of:

- `:NOT_POSITIVE_DEFINITE` if the left-hand-side matrix `A` is found to be not
  positive definite;

- `:TOO_MANY_ITERATIONS` if the maximum number of iterations have been reached;

- `:F_TEST_SATISFIED` if convergence occured because the function reduction
  satisfies the criterion specified by `ftol`;

- `:G_TEST_SATISFIED` if convergence occured because the gradient norm
  satisfies the criterion specified by `gtol`;

- `:X_TEST_SATISFIED` if convergence occured because the norm of the variation
  of variables satisfies the criterion specified by `xtol`.

Method [`ConjugateGradient.reason`](@ref) may be called to get a textual
explanation about the returned status.

"""
function solve!(x::V, A, b::V, ctx::Context{V},
                io::IO = devnull) where {V}
    # Run the unpreconditioned version of the algorithm.
    return solve!(x, A, b, Identity(), ctx, io)
end

function solve!(x::V, A, b::V, M, ctx::Context{V},
                io::IO = devnull) where {V}
    # Get workspace variables.
    precond = !isa(M, Identity)
    verbose = !isa(io, Base.DevNull)
    p, q, r = ctx.p, ctx.q, ctx.r
    z = (precond ? ctx.z : ctx.r)
    if precond && z === r
        error("workspace variables Z must be different from R with ",
              "a preconditioner")
    end

    # Enforce types of some variables (FIXME: This will not work for BigFloat).
    local rho::Float64, oldrho::Float64
    local alpha::Float64, gamma::Float64
    local psi::Float64, psimax::Float64
    local gtest::Float64

    # Initialize local variables.
    rho = 0.0
    psi = 0.0
    psimax = 0.0
    xtest = (ctx.xatol > 0 || ctx.xrtol > 0)
    if verbose
        t0 = time()
    end

    # Conjugate gradient iterations.
    k = 0
    while true
        # Is this the initial or a restarted iteration?
        restart = (k == 0 || ctx.restart > 0 && rem(k, ctx.restart) == 0)

        # Compute residuals and their squared norm.
        if restart
            # Compute residuals.
            if k > 0 || vnorm2(x) != 0
                # Compute r = b - A*x.
                vcombine!(r, 1, b, -1, vmul!(r, A, x))
            else
                # Spare applying A since x = 0.
                vcopy!(r, b)
            end
        else
            # Update residuals.
            vupdate!(r, -alpha, q)
        end
        if precond
            # Apply preconditioner.
            vmul!(z, M, r) # z = M*r
        end
        oldrho = rho
        rho = vdot(r, z) # rho = ‖r‖_M^2
        if k == 0
            gtest = tolerance(ctx.gatol, ctx.gatol, sqrt(rho))
        end
        if verbose
            t = (time() - t0)*1E3 # elapsed time in ms
            if precond
                if k == 0
                    print(io,
                          "# Iter.   Time (ms)     Δf(x)       ‖∇f(x)‖  ",
                          "   ‖∇f(x)‖_M\n",
                          "# --------------------------------------------",
                          "-------------\n")
                end
                @printf(io, "%7d %11.3f %12.4e %12.4e %12.4e\n",
                        k, t, psi, vnorm2(r), sqrt(rho))
            else
                if k == 0
                    print(io,
                          "# Iter.   Time (ms)     Δf(x)       ‖∇f(x)‖\n",
                          "# --------------------------------------------\n")
                end
                @printf(io, "%7d %11.3f %12.4e %12.4e\n",
                        k, t, psi, sqrt(rho))
            end
        end
        if sqrt(rho) ≤ gtest
            # Normal convergence in the gradient norm.
            verbose && println(io, "# Convergence in the gradient norm.")
            return :G_TEST_SATISFIED
        end
        if k ≥ ctx.maxiter
            verbose && println(io, "# Too many iteration(s).")
            return :TOO_MANY_ITERATIONS
        end

        # Compute search direction.
        if restart
            # Restarting or first iteration.
            vcopy!(p, z)
        else
            # Apply recurrence.
            beta = rho/oldrho
            vcombine!(p, 1, z, beta, p)
        end

        # Compute optimal step size.
        vmul!(q, A, p)
        gamma = vdot(p, q)
        if !(gamma > 0)
            verbose && println(io, "# Operator is not positive definite.")
            return :NOT_POSITIVE_DEFINITE
        end
        alpha = rho/gamma

        # Update variables and check for convergence.
        vupdate!(x, +alpha, p)
        psi = alpha*rho/2  # psi = f(x_{k}) - f(x_{k+1}) ≥ 0
        psimax = max(psi, psimax)
        if psi ≤ tolerance(ctx.fatol, ctx.frtol, psimax)
            # Normal convergence in the function reduction.
            verbose && println(io, "# Convergence in the function reduction.")
            return :F_TEST_SATISFIED
        end
        if xtest && alpha*vnorm2(p) ≤ tolerance(ctx.xatol, ctx.xrtol, x)
            # Normal convergence in the variables.
            verbose && println(io, "# Convergence in the variables.")
            return :X_TEST_SATISFIED
        end

        # Increment iteration number.
        k += 1
    end
end

"""

Given absolute and relative tolerances `atol` and `rtol` (both finite and
nonnegative), the calls:

    tolerance(atol, rtol, val) -> max(0, atol, rtol*abs(val))
    tolerance(atol, rtol, arr) -> max(0, atol, rtol*vnorm2(arr))

yield the tolerances for the scalar `val` or for the array `arr`. The result is
nonnegative.

If `rtol ≤ 0`, the computation of `vnorm2(arr)` is not performed.

"""
tolerance(atol::Real, rtol::Real, val::Real) =
    tolerance(promote(atol, rtol, val)...)

tolerance(atol::T, rtol::T, val::T) where {T<:Real} =
    max(zero(T), atol, rtol*abs(val))

function tolerance(atol::Real, rtol::Real,
                   arr::AbstractArray{T}) where {T<:AbstractFloat}
    # NOTE: This is to ensure type stability.
    val = (rtol > 0 ? vnorm2(arr)::T : zero(T))
    tolerance(atol, rtol, val)
end

"""
    bad_argument(args...)

throws an `ArgumentError` exception with error message given by `args...`
converted into a string.

"""
bad_argument(msg::AbstractString) = throw(ArgumentError(msg))
@noinline bad_argument(args...) = bad_argument(string(args...))

end # module
