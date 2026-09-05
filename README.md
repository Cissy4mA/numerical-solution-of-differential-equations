# Numerical Solution of Differential Equations — MATLAB

> Classical finite-difference and linear-multistep methods for ODE initial-value problems and second-order PDEs (parabolic, elliptic, hyperbolic), with step-size and scheme comparisons, error analysis, and convergence observations.

## Overview

| Category | Prototype | Numerical methods | Directory |
|----------|-----------|-------------------|-----------|
| ODE IVP | `y' = f(t, y), y(t0) = y0` | 2nd-order explicit Adams, trapezoidal (implicit Adams), Euler; step-size and scheme error comparison | `ODE/` |
| Parabolic PDE | `u_t = a u_xx + f(x,t)` | Forward (explicit), Crank–Nicolson, upwind; Thomas/chasing for tridiagonal systems | `Parabolic/` |
| Elliptic PDE | Poisson `u_xx + u_yy = f(x,y)` | Five-point stencil with Kronecker-product block matrix; grids `h = 1/64` and `1/128` | `Elliptic/` |
| Hyperbolic PDE | Wave equation `u_tt = a^2 u_xx` | Explicit, implicit, and upwind schemes; block-tridiagonal solver | `Hyperbolic/` |

## Directory Structure

```
numerical-solution-of-differential-equations/
├── ODE/            # ODE initial-value problems
├── Parabolic/      # Parabolic PDEs
├── Elliptic/       # Elliptic PDEs
├── Hyperbolic/     # Hyperbolic PDEs
├── README.md
└── LICENSE
```

## 1. ODE Initial-Value Problems (`ODE/`)

**Problem.** Solve the nonlinear IVP

```math
y'(t) = 4 t \sqrt{y(t)}, \quad y(0) = 1, \quad t \in [0, 2]
```

whose exact solution is `y(t) = (t^2 + 1)^2`.

**Methods.**

- **2nd-order explicit Adams (Bashforth).** Two-step predictor using previously computed slopes:

```math
y_{n+1} = y_n + \frac{h}{2}\bigl(3 f(t_n, y_n) - f(t_{n-1}, y_{n-1})\bigr)
```

- **Trapezoidal (implicit).** One-step implicit Adams–Moulton rule, solved algebraically at each step:

```math
y_{n+1} = y_n + \frac{h}{2}\bigl(f(t_n, y_n) + f(t_{n+1}, y_{n+1})\bigr)
```

- **Euler.** Used only to generate the starting value `y_1` for the Adams scheme.

**Findings.**

- Halving the step size from `h = 0.1` to `h = 0.05` visibly improves the Adams trajectory, confirming the expected convergence order.
- At the same step size `h = 0.1`, the trapezoidal result is almost indistinguishable from the exact solution, while the explicit Adams still deviates.
- Cost trade-off: the trapezoidal method requires solving an implicit equation per step (higher CPU time), whereas explicit Adams is cheaper per step but less accurate for the same `h`.

**Files.**

- `Adams_2.m` — explicit Adams with symbolic exact-solution overlay
- `h1h2Adams.m` — step-size comparison (`h = 0.1` vs `0.05`) for explicit Adams
- `ladder_adams.m` — trapezoidal vs explicit Adams at the same step size

## 2. Parabolic PDEs (`Parabolic/`)

**Problem.** 1-D heat-type equation with source term on `(x, t) \in (0, 1) \times (0, T]`:

```math
u_t = a u_{xx} + f(x, t)
```

with given initial and Dirichlet boundary data. The example implemented uses exact solution `u(x,t) = e^{-2\pi^2 t} \cos(\pi x) + 1 - \cos(t)`.

**Methods.**

- **Forward (explicit) scheme.** First-order in time, second-order in space; conditionally stable under the mesh-ratio restriction `r = a \tau / h^2 \le 1/2`.
- **Crank–Nicolson.** Second-order in both space and time, unconditionally stable; each time step requires solving a tridiagonal system, handled by the Thomas (chasing) algorithm.
- **Upwind.** A one-sided spatial discretization for convection-dominated or asymmetric settings.

**Implementation notes.** The coefficient matrix for a `(m-1)`-node spatial discretization is built from block-tridiagonal Kronecker products of the identity and the standard second-difference matrix. The right-hand-side vector incorporates the source term plus the initial layer.

**Findings.**

- Forward difference is easy to implement but is constrained by a strict stability limit on the time step.
- Crank–Nicolson is more robust: larger time steps are allowed while preserving second-order accuracy, making it preferable for parabolic problems where efficiency matters.

**Files.**

- `GrankNicolsonParabolic.m` — Crank–Nicolson driver (calls `Chasing3`)
- `Chasing3.m` — Thomas algorithm for tridiagonal systems
- `paowu.m` — forward (explicit) scheme
- `paowu_GN.m` — Crank–Nicolson scheme
- `paowu_xiangqian.m` — upwind scheme
- `testParabolic.m` — integrated test and animation frame capture

## 3. Elliptic PDEs (`Elliptic/`)

**Problem.** 2-D Poisson equation with zero Dirichlet boundary conditions on the unit square:

```math
-\Delta u = f(x, y), \quad (x,y) \in (0,1)^2
```

with the test case exact solution `u(x,y) = e^{\pi(x+y)} \sin(\pi x) \sin(\pi y)`.

**Method.** Standard five-point finite-difference stencil:

```math
\frac{4 u_{i,j} - u_{i+1,j} - u_{i-1,j} - u_{i,j+1} - u_{i,j-1}}{h^2} = f_{i,j}
```

The resulting linear system is assembled as a block-tridiagonal matrix `DCD` using Kronecker products: `main_diag = I \otimes C`, `off_diagonals = E \otimes D`, where `C` is the 1-D second-difference matrix and `D` is the coupling matrix for the `y` direction. The interior nodal values are obtained by direct inversion and then padded back onto the full grid.

**Findings.**

- At `h = 1/64` and `h = 1/128` the numerical solution tracks the exact solution closely across all interior nodes.
- Finer grids reduce the truncation error as expected from the `O(h^2)` spatial discretization.
- Visualizing the solution surface (`mesh` plots) confirms the boundary conditions and the qualitative shape of the exact solution.

**Files.**

- `tuoyuan.m` — main five-point Poisson solver
- `tuoyuan_5point.m` — five-point stencil with `h = 1/64` and `1/128` comparisons, error plots, and 3-D surface visualization

## 4. Hyperbolic PDEs (`Hyperbolic/`)

**Problem.** 1-D wave equation (mixed initial-boundary-value problem):

```math
u_{tt} = a^2 u_{xx}, \quad x \in (0,1), \; t \in (0,2]
```

with initial displacement `u(x,0) = 2 \sin(\pi x)` and exact solution `u(x,t) = \sin\pi(x-t) + \sin\pi(x+t)`.

**Methods.**

- **Explicit scheme.** Three-time-level stencil:

```math
u_j^{k+1} = 2(1-r^2) u_j^k + r^2\bigl(u_{j+1}^k + u_{j-1}^k\bigr) - u_j^{k-1}, \quad r = a \tau / h
```

Stable when `r \le 1` (CFL condition).

- **Implicit scheme.** The future time level `u^{k+1}` is coupled into the stencil, producing a block-tridiagonal linear system solved by the chasing method. This removes the CFL restriction but increases per-step cost.

- **Upwind scheme.** One-sided spatial discretization for transport-dominated behavior.

**Findings.**

- Both explicit and implicit schemes produce numerical solutions close to the exact traveling-wave solution.
- At fixed sampling times, the explicit and implicit results are visually similar; the explicit scheme is cheaper per step but requires satisfying the CFL condition.

**Files.**

- `shuangqu.m` / `shuangqu_xian.m` — explicit and implicit wave-equation solvers
- `Chasing.m` / `Chasing3.m` — tridiagonal and block-tridiagonal chasing solvers
- `Hyperbolic2.m` — additional format implementation
- `testHyperbolic.m` / `testHyperbolic2.m` — integrated tests and comparison plots

## Running the Code

- Requires **MATLAB** with basic numerical and plotting support; no external toolboxes are needed.
- Each `.m` file is self-contained. Scripts that start with `test*.m` are intended as entry points that run a complete method and produce comparison plots.
- For the Parabolic and Hyperbolic implicit schemes, make sure the corresponding `Chasing*.m` files are on the MATLAB path.

## License

MIT — see [LICENSE](LICENSE).
