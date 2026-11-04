# Equations and numerical methods

This self-contained specification describes the equations, initial data,
parameters and numerical diagnostics implemented by this software. It is the
reference for reproducing the supplied synthetic experiments.

## Model

Positions are ordered real lifts of points on a circle of length one:

$$x_{i+N}=x_i+1,\qquad g_i=N(x_{i+1}-x_i),\qquad N^{-1}\sum_i g_i=1.$$

The particle equation is

$$\dot x_i=U^N(x_i)+\nu_N N\{A(g_i)-A(g_{i-1})\},$$

where

$$M^N(x)=1+N^{-1}\sum_j K(x-x_j)\Phi(g_j),\quad
F^N(x)=N^{-1}\sum_j W'(x-x_j),\quad U^N=-M^NF^N.$$

The coefficients are

$$K(x)=1+\tfrac12\cos(2\pi x),\quad W(x)=-\cos(2\pi x),\quad
\Phi(g)=e^{-g},\quad A(g)=g+\beta g^2.$$

The main experiment has `beta = 1/2`. The affine comparison has `beta = 0`.
Sine and cosine moments evaluate the nonlocal fields in linear work per state
evaluation; the fast checks compare them with direct pairwise sums.

The arranged-cell reference satisfies

$$X_t=U(X),\qquad G_t=U_x(X)G+\lambda\partial_{\theta\theta}A(G),\qquad
q=\int_0^1G\,d\theta=X_a,$$

with

$$H(a)=\int_0^1e^{-G(a,\theta)}\,d\theta,\quad
M(x)=1+\int_0^1K(x-X(a))H(a)\,da,\quad
F(x)=\int_0^1W'(x-X(a))\,da,\quad U=-MF.$$

## Initial preparation and settings

The density is `rho0(x) = 1 + cos(2*pi*x)/3`. Its quantile map solves
`X0 + sin(2*pi*X0)/(6*pi) = a`, and `q0 = 1/rho0(X0)`.
The two cells are

$$G_0^{(k)}(a,\theta)=q_0(a)[1+\tfrac14\cos(2\pi k\theta)],\qquad k=1,2.$$

They have identical value distributions at each macroscopic label. Their cell
derivatives differ, so positive cell diffusion can distinguish their evolution.

Particles sample `q0(i/N) * [1 + cos(2*pi*k*i/p)/4]`; these gaps are normalized
by their discrete mean. Positions start at zero and are recovered by cumulative
gap sums divided by `N`. The two finite preparations need not have identical
initial empirical position measures; the zero-diffusion control measures this
finite-size effect.

| Setting | Value |
| --- | --- |
| `n`, `N`, `p` | `n = 8,16,32`; `N = n^3`; `p = n` |
| Particle diffusivity | `nu_N = lambda * (p/N)^2` |
| Main/control diffusivity | `lambda = 1 / 0` |
| Final time | `T = 0.02` |
| Uniform output mesh | 101 times, including zero and `T` |
| Additional response times | `1e-3, 5e-4, 2e-4, 1e-4, 5e-5` |
| Fine cell grid | 512 macroscopic labels × 64 fast labels |
| Coarse cell grid | 256 × 32 |
| Fine/coarse cell RK4 maximum step | `5e-6 / 2e-5` |
| Particle RK4 maximum step | `1e-5`, refined to `5e-6` |
| Affine main comparison | 1024 macroscopic labels |
| Continuum transport quadrature | 8192 labels, refined to 16384 |

RK4 steps are shortened to land exactly on requested output times. The cell
diffusion uses Fourier collocation in the fast label, with pointwise evaluation
of `A(G)`; there is no dealiasing filter. Label integrals use periodic trapezoidal
quadrature. Spatial refinement, an independent finite-difference discretization,
and exact affine solutions check this discretization for the stated smooth data.
The hard-coded steps and grids are calibrated for these experiments, not a
general parameter-sweep solver. Changing diffusivity, cell frequency or grid
size requires reassessing stability and accuracy.

## Reported quantities

The joint error at each output time is `E_N + D_N`, where `E_N` is the maximum
absolute gap error and `D_N` the maximum difference of lifted positions. Periodic
cubic interpolation samples the cell at the actual labels and fast phases;
positions are interpolated through their periodic displacement `X-a`.
For the `k=1` preparation, the reported error is the maximum of `E_N + D_N`
over the 101 uniformly spaced output times.
Along `N=n^3`, both `p/N` and `p^(-2)` equal `N^(-2/3)`, giving the plotted guide.

The observable is

$$J_k(t)=\int_0^1\psi(X^{(k)}(t,a))\,da,\qquad
\psi(x)=-\tfrac16\cos(2\pi x).$$

Its difference is evaluated through a sine-product identity to reduce numerical
cancellation. Independent initial-data quadrature predicts
`(J2-J1)/t^2 -> C2/2`: approximately 0.741984411175 for nonlinear diffusion and
0.365253559927 for affine diffusion. The zero-time coefficients computed from
this integral are stored in `theoretical_coefficients.csv`. The response data
cover positive times.

### Explicit response coefficient

Let `G0 = G0^(1)` be the frequency-one cell and define

$$F_0(x)=\frac{\pi}{3}\sin(2\pi x),\qquad
D_\beta(a)=\int_0^1 e^{-G_0(a,\theta)}[1+2\beta G_0(a,\theta)]
                 |\partial_\theta G_0(a,\theta)|^2\,d\theta.$$

Here `partial_theta G0 = -2*pi*delta*q0(a)*sin(2*pi*theta)`, with `delta=1/4`.
The coefficient computed by `simulation.coefficient(beta)` is

$$C_2=3\int_0^1 F_0(X_0(a))^2
       \left[\int_0^1 K(X_0(a)-X_0(z))D_\beta(z)\,dz\right]da.$$

The factor 3 is `k^2-1` for `k=2`; the implemented calculation takes `lambda=1`
and unit mobility coupling. For general fixed `lambda`, this expression is
multiplied by `lambda`. The inner convolution is evaluated by zeroth and first
sine/cosine moments, and both integrations use periodic quadrature. The
independent implementation changes to physical coordinates with `da=rho0(x) dx`
and uses adaptive quadrature.

The formula follows directly from the model. Differentiating
`H=integral Phi(G) dtheta` and integrating the diffusion term by parts gives

$$H_t=b\int_0^1G\Phi'(G)\,d\theta
       -\lambda\int_0^1\Phi''(G)A'(G)|G_\theta|^2\,d\theta.$$

At time zero the two frequencies have the same value law, map and velocity.
The first term is therefore the same, while the second cell's gradient term is
four times the first one's. Thus `H_t^(2)-H_t^(1)=-3*lambda*D_beta` initially.
The force time derivatives agree initially, so the difference of the map
accelerations is `3*lambda*F0*(K*D_beta)`. In `J''`, the velocity-squared terms
cancel and `psi'=F0`; integration gives `J2''(0)-J1''(0)=C2`. Since the initial
observable difference and its first derivative vanish, Taylor expansion yields
`J2-J1 = (C2/2)*t^2 + O(t^3)`. This explains both the coefficient and the
factor one half in the tables without a separate publication.

Circular `W1` uses the geodesic distance, not the distance on a line and not a
fixed label matching. For equal-weight empirical measures, allowing different
numbers of atoms, the implementation minimizes the integral of
`abs(F1-F2-c)` over a constant `c`, using a length-weighted median of the CDF
difference. This is exact for the supplied atoms up to floating-point arithmetic.
Continuum `W1` additionally incurs sampling error: smooth periodic displacement
maps are resampled before the empirical distance is evaluated. Linear-program
transport solutions independently check the atomic formula, including repeated
atoms, circle wrap-around, unequal counts, rotations and symmetry.

## Independent audit

`validation.py` uses conservative gap variables plus the mean lifted position,
complex-exponential field formulas, full complex FFTs, and SciPy's DOP853
integrator. This checks the production position formulation, trigonometric
moments, real FFTs, and RK4 integration through different implementations.
The independence is numerical and algebraic; both implementations target the
same model and share the prescribed observation times and transport routine.
The transport routine is checked separately against linear programming.

The DOP853 tolerances are `rtol=2e-12`, `atol=2e-14`. Diffusion-scale maximum steps
are `1/(4*lambda*p^2)` for particles and `1/[lambda*(pi*m)^2]` for `m` fast cell
points; at zero diffusivity those caps are omitted. The initial step is at most
`1e-6`. These caps control high-frequency modes, including dense-output errors.
Do not remove them on the assumption that tolerance control alone guarantees
stability. A 100-fold tolerance change is also checked; it may have little effect
when the diffusion step cap is active.

Additional checks include separate label/cell refinements, centered finite
differences at `m=16,32,64,128`, exact affine heat evolution, gap positivity,
mass, `q=X_a`, the cell-law and variance identities, and the variance decay bound.
Physical-coordinate adaptive quadrature independently computes the response
coefficient. Small-time integration at times down to `5e-6` tests its limit.

Fast relaxation uses `lambda = 1,2,4,8,16,32`, a 64 × 32 cell grid, and an
independently integrated reduced model with `H=exp(-q)` and `q_t=U_x(X)q`.
The reported error is the time-mesh maximum of `max_a|X-Xbar| + max_a|q-qbar|`.
A 64-point fast-label run checks the `lambda=32` result. For
`r=G/q`, the variance `integral (r-1)^2 dtheta` is compared with
`(1/32)*exp(-8*pi^2*lambda*t)`.

`affine.py` reproduces the supplementary affine experiment on `T=0.01` for
`N=64,512,4096`, using the exact heat-cell form and a Bessel-series coefficient
check. `certificate.py` evaluates rational bounds with `fractions.Fraction` and
verifies the separate affine observation time `729/(8192*10^9)` exactly.

### What the rational certificate checks

The certificate concerns only the smooth positive affine model (`beta=0`,
`lambda=1`, `k=1,2`) on `0 <= t <= T0=1e-3`. It is a conservative bound, not
a fitted lifetime of the observed separation. All norms below are suprema
over both solutions, their labels and this interval.

Let `S_j` bound `|partial_theta^j G|`, `K_j` bound the spatial derivatives of
`K`, and `w_j` bound those of `W`. Because `Phi=exp(-g)` and `g>0`, bounds for
`|Phi|`, `|Phi'|` and `|Phi''|` are all one. The field bounds are
`|U_x| <= B=12*pi^2` and `|U_xx| <= B2=32*pi^3`. From the affine heat formula
and `q_t=U_x(X)q`,

$$S_0\leq\frac{15}{8}e^{BT_0},\qquad
S_j\leq\frac38 e^{BT_0}(4\pi)^j\quad(1\leq j\leq4).$$

The code uses rational upper bounds by replacing `pi` with `r=22/7`:

$$K_0=\frac32,\ K_1=r,\ K_2=2r^2,\quad w_j=(2r)^j,
\quad B=12r^2,\ B_2=32r^3,\quad S_0=\frac{15}{7},\qquad
S_j=\frac37(4r)^j.$$

These `S_j` values follow from `B*T0 < 3/25` and
`exp(3/25) <= 1/(1-3/25) = 25/22 < 8/7`.
The remaining variables in `certificate.py` have these meanings:

| Code variable | Quantity bounded |
| --- | --- |
| `M0`, `M1` | `M`, `M_x` |
| `V0`, `V2`, `V3` | `X_t`, `X_tt`, `X_ttt` |
| `Q1`, `Q11`, `Q12`, `Q2` | `G_t`, its first/second cell derivatives, and `G_tt` |
| `R0`, `R1`, `Rtt` | `M_t`, `M_xt`, `M_tt` |
| `T0`, `T1`, `T2` in the code | `U_t`, `U_xt`, `U_tt` |

The code variable `T0` in this table is a derivative bound, not the interval
endpoint `T0=1e-3` above. Likewise, `V0` here is a speed bound, not cell variance.
For affine diffusion the recurrence simplifies to

$$\begin{aligned}
M_0&=1+K_0,& M_1&=K_1,& V_0&=M_0w_1,\\
Q_1&=BS_0+S_2,& Q_{1,1}&=BS_1+S_3,& Q_{1,2}&=BS_2+S_4,\\
R_0&=K_1V_0+K_0Q_1,& R_1&=K_2V_0+K_1Q_1,\\
T_0&=R_0w_1+M_0w_2V_0,\\
T_1&=R_1w_1+R_0w_2+M_1w_2V_0+M_0w_3V_0,\\
Q_2&=(T_1+B_2V_0)S_0+BQ_1+Q_{1,2},& V_2&=T_0+BV_0,\\
R_{tt}&=K_2V_0^2+K_1V_2+2K_1V_0Q_1+K_0(Q_1^2+Q_2),\\
T_2&=R_{tt}w_1+2R_0w_2V_0+M_0(w_3V_0^2+w_2V_2),\\
V_3&=T_2+2T_1V_0+B_2V_0^2+BV_2.
\end{aligned}$$

These inequalities come from differentiating `U=-MF`, the source integrals,
the affine cell equation and `X_t=U(X)`, then applying the triangle inequality.
For example `X_tt=U_t+U_x*X_t` gives `V2`, and
`X_ttt=U_tt+2*U_xt*X_t+U_xx*X_t^2+U_x*X_tt` gives `V3`.
Because `||psi^(j)|| <= w_j`, the third derivative of the observable difference
is bounded by

$$L_3=2(w_3V_0^3+3w_2V_0V_2+w_1V_3)<10^9.$$

This last inequality is evaluated with exact `fractions.Fraction` arithmetic.
The initial coefficient has the independent lower bound
`C2 >= (3*pi^4/512)*exp(-15/8) > 243/4096`, using `pi>3` and `exp(15/8)<8`.
Taylor's inequality consequently gives

$$J_2(t)-J_1(t)\geq\frac{C_2}{2}t^2-\frac{L_3}{6}t^3
\geq\frac{C_2}{4}t^2\quad\text{for }0\leq t\leq
\frac{729}{8192\cdot10^9}.$$

The certified time is about `8.899e-11`, lies inside `T0=1e-3`, and is distinct
from the simulations' `T=0.02`. The certificate checks this explicit bound; it
does not certify the numerical trajectories at the longer simulation times.

## Interpretation

These checks support the predicted behavior for the stated examples. They do
not prove theorems, certify continuum suprema, establish uniformity in untested
parameter limits, or validate a traffic, pedestrian, cellular or granular model.
The tiny rationally certified time and the computational horizon `T=0.02` have
different purposes and must not be conflated.
