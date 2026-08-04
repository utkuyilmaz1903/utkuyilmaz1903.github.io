+++
title = "What Non-Uniform WENO Actually Buys Us"
description = "Quantitative showcase of non-uniform WENO-5 on viscous shocks and multi-domain front compression in MethodOfLines.jl."
tags = ["gsoc", "julia", "sciml", "pde", "weno"]
hasmath = true
hascode = true
date = Date(2026, 8, 1)
rss_title = "What Non-Uniform WENO Actually Buys Us"
rss_description = "Two physical scenarios showing where clustered and multi-domain non-uniform WENO pays off against exact solutions."
+++

@@blog-post

# What Non-Uniform WENO Actually Buys Us

~~~
<header class="article-header">
  <p class="article-series">GSoC 2026 · MethodOfLines.jl · Final DevLog</p>
  <p class="article-meta">August 2026 · Quantitative showcase</p>
  <div class="tag-row">
    <span class="tag-pill">WENO-5</span>
    <span class="tag-pill">Non-uniform grids</span>
    <span class="tag-pill">SciML</span>
  </div>
</header>
~~~

~~~
<p class="lead">Earlier posts in this series derived the mathematics and built the implementation. This final post answers a practical question: <em>what does that infrastructure buy in a physical simulation?</em></p>
~~~

~~~
<div class="summary-box">
  <p class="summary-title">At a glance</p>
  <ul>
    <li><strong>Scenario A</strong> — A viscous shock layer on $[-1,1]$: clustered WENO at $N=129$ beats uniform WENO at $N=257$ by roughly two orders of magnitude.</li>
    <li><strong>Scenario B</strong> — A compressing interfacial front across two mismatched non-uniform subdomains: machine-precision seam continuity against an exact characteristic solution.</li>
    <li><strong>Takeaway</strong> — Under-resolving a steep feature on a uniform grid is a failure regime, not a convergence problem. Grid design is as important as the operator.</li>
  </ul>
</div>
~~~

## Introduction

My GSoC 2026 project for [MethodOfLines.jl](https://github.com/SciML/MethodOfLines.jl) is complete. The goal was to bring native non-uniform grid support to the WENO-5 advection scheme. Earlier posts covered the [smoothness indicators](/weno-smoothness/) and the [boundary architecture](/nonuniform-weno-boundaries/).

Many PDE problems contain localized steep gradients: shock-like transitions, boundary layers, or interfacial fronts. On a uniform grid, resolving one of these features means placing enough points everywhere to resolve it locally. Most of the grid then sits in regions where the solution is nearly flat. The cost scales with the global mesh, not with the feature.

There is a second issue. Under-resolving a steep gradient on a uniform grid is not a small-error regime. When the layer width falls below the grid scale, the semi-discretization can become unstable. The simulation does not converge slowly; it fails. Non-uniform WENO addresses both problems by concentrating resolution where the physics demands it, while keeping a high-order, non-oscillatory reconstruction on stretched meshes.

The two scenarios below quantify that difference against exact solutions.

<hr class="section-rule">

@@scenario-block
<span class="scenario-label">Scenario A</span>

## The Stationary Viscous Shock

The first problem is the viscous Burgers equation on $[-1, 1]$ with viscosity $\nu = 2 \times 10^{-3}$:

$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2},
\qquad u(t, \mp 1) = \pm 1.
$$

Nonlinear advection steepens the profile; diffusion regularizes it. The equation admits the exact steady state

$$
u_\infty(x) = -\tanh\left(\frac{x}{2\nu}\right),
$$

an interior layer of width $\mathcal{O}(\nu)$ centered at $x = 0$. With $\nu = 2 \times 10^{-3}$, the layer is roughly $0.01$ wide. It is a standard finite-difference prototype for a steep interfacial gradient at a known, fixed location.

The test protocol is deliberate. We initialize from the exact steady profile and integrate to $t = 1$. A capable scheme holds the layer in place with small error and no spurious oscillation. A poor scheme smears it, shifts it, or becomes numerically unstable.

**We compare three configurations:**

~~~
<div class="config-grid">
  <div class="config-card">
    <span class="config-name">First-order upwind</span>
    <span class="config-spec">Uniform · $N = 129$</span>
  </div>
  <div class="config-card">
    <span class="config-name">WENO-5</span>
    <span class="config-spec">Uniform · $N = 257$</span>
  </div>
  <div class="config-card config-card--highlight">
    <span class="config-name">WENO-5</span>
    <span class="config-spec">Clustered · $N = 129$</span>
  </div>
</div>
~~~

The clustered grid is built from a resolution density peaked at the layer. Points are placed at uniform quantiles of that density, so cell widths vary smoothly across the domain. In MethodOfLines.jl, the full setup is a standard `PDESystem` discretized with `WENOScheme()` on a custom grid vector:

~~~
<p class="code-caption">Listing 1 — Non-uniform WENO discretization via a density-clustered grid</p>
~~~

```julia
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets
using OrdinaryDiffEqSSPRK: SSPRK33

@parameters t x
@variables u(..)
Dt = Differential(t)
Dx = Differential(x)

ν = 2.0e-3
u_steady(x) = -tanh(x / (2ν))

function grid_from_density(a, b, n, ρ; nsamples = 5001)
    xs = range(a, b, length = nsamples)
    cdf = cumsum(ρ.(xs))
    cdf = (cdf .- cdf[1]) ./ (cdf[end] - cdf[1])
    xg = map(range(0, 1, length = n)) do l
        k = searchsortedfirst(cdf, l)
        k <= 1 && return float(xs[1])
        θ = (l - cdf[k - 1]) / (cdf[k] - cdf[k - 1])
        return xs[k - 1] + θ * (xs[k] - xs[k - 1])
    end
    xg[1] = a
    xg[end] = b
    return xg
end

cluster_density(x) = 1 + 30 * exp(-x^2 / (2 * 0.02^2))
xg = grid_from_density(-1.0, 1.0, 129, cluster_density)

eq = Dt(u(t, x)) ~ -u(t, x) * Dx(u(t, x)) + ν * Dx(Dx(u(t, x)))
bcs = [u(0.0, x) ~ u_steady(x), u(t, -1.0) ~ 1.0, u(t, 1.0) ~ -1.0]
domains = [t ∈ Interval(0.0, 1.0), x ∈ Interval(-1.0, 1.0)]
@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

disc = MOLFiniteDifference([x => xg], t; advection_scheme = WENOScheme())
prob = discretize(pdesys, disc)
```

The grid vector `xg` is the only non-standard input. Everything else follows the usual MethodOfLines workflow.

@@

~~~
<figure class="blog-figure">
  <img src="/assets/01_layer_after_t1.png" alt="Profile zoom after t = 1">
  <figcaption>
    <span class="fig-label">Figure 1</span>
    Profile zoom after $t = 1$. Clustered WENO with $N = 129$ tracks the exact steady state. Uniform WENO at $N = 257$ remains non-oscillatory but visibly inaccurate. First-order upwind is dominated by numerical viscosity and resolves the wrong layer width.
  </figcaption>
</figure>
~~~

~~~
<p class="table-caption">Table 1 — Relative $L^2$ error after $t = 1$</p>
~~~

~~~
<div class="table-wrapper">
<table class="results-table">
  <thead>
    <tr>
      <th>Configuration</th>
      <th>Relative $L^2$ error</th>
      <th>Notes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Upwind-1, uniform $N = 129$</td>
      <td>$\approx 1.6 \times 10^{-2}$</td>
      <td>Numerical viscosity $\gg$ physical $\nu$; no convergence in this range</td>
    </tr>
    <tr class="row-unstable">
      <td>WENO, uniform $N = 129$</td>
      <td>—</td>
      <td>Unstable: sub-cell layer aborts the integration</td>
    </tr>
    <tr>
      <td>WENO, uniform $N = 257$</td>
      <td>$\approx 4.6 \times 10^{-2}$</td>
      <td>Stable, but ~2 orders of magnitude worse than clustered</td>
    </tr>
    <tr class="row-best">
      <td>WENO, clustered $N = 129$</td>
      <td>$\approx 1.3 \times 10^{-4}$</td>
      <td>Overshoot at the $10^{-8}$ level</td>
    </tr>
  </tbody>
</table>
</div>
~~~

~~~
<div class="insight-box">
  <p class="insight-label">Key result</p>
  <p>The clustered grid with <strong>half the degrees of freedom</strong> of the stable uniform WENO run delivers roughly two to three orders of magnitude lower error. The operator is the same; the grid design is not.</p>
</div>
~~~

~~~
<figure class="blog-figure">
  <img src="/assets/02_accuracy_per_grid_point.png" alt="Relative L2 error versus grid points">
  <figcaption>
    <span class="fig-label">Figure 2</span>
    Relative $L^2$ error versus grid points on a log-log scale. Clustered WENO converges cleanly across a self-similar refinement family. Uniform WENO appears as a single point at $N = 257$ because the $N = 129$ run is unstable. Upwind does not converge until $N \gtrsim 1000$, when its numerical viscosity finally drops below the physical $\nu$.
  </figcaption>
</figure>
~~~

<hr class="section-rule">

@@scenario-block
<span class="scenario-label">Scenario B</span>

## Interfacial Front Compression

The second problem moves from a single domain to a two-domain setup inspired by transport in a flow battery. Electrolyte enters from an open channel and flows into a porous electrode. Interstitial velocity decreases as porosity drops, and any concentration front carried by the flow compresses as it crosses the transition.

We model along-flow advection with a spatially varying slowness field

$$
s(x) = 1 + \left(\frac{1}{v_2} - 1\right) \frac{1 + \tanh((x - x_v)/\delta_v)}{2},
$$

which ramps smoothly from $1$ to $1/v_2 = 2$ across a transition zone. The exact solution follows from characteristics:

$$
c(x, t) = S(t - \tau(x)), \qquad \tau(x) = \int_0^x s(\xi)\, d\xi,
$$

where $S$ is a smooth inlet signal. A front of temporal width $w_t$ has spatial width $v(x)\, w_t$. It enters at width $0.02$ and leaves the transition at width $0.01$.

The domain is split at $x = 1/2$ into two subdomains with independently generated non-uniform grids. Domain 2 uses roughly twice the point density, because the compressed front lives there. Continuity at the seam is enforced algebraically. The interface is computational: the underlying velocity field $v(x) = 1/s(x)$ is continuous across the full domain.

@@

~~~
<figure class="blog-figure">
  <img src="/assets/03_front_compression.png" alt="Front compression across the porosity transition">
  <figcaption>
    <span class="fig-label">Figure 3</span>
    A steep front arrives at the seam at $t \approx 0.67$, crosses onto a finer grid without a visible kink, and emerges at $t = 1.3$ twice as steep. Numerical points track the exact characteristic solution on both sides.
  </figcaption>
</figure>
~~~

**Against the exact solution:**

~~~
<ul class="results-checklist">
  <li>
    <span class="check-title">Seam continuity</span>
    <span class="check-detail">Machine precision. The interface identification is algebraic, so $|c_1 - c_2|$ at the seam is zero to floating-point tolerance.</span>
  </li>
  <li>
    <span class="check-title">Boundedness</span>
    <span class="check-detail">The solution remains in $[0, 1]$ to within $10^{-3}$ while the steep front crosses grids.</span>
  </li>
  <li>
    <span class="check-title">Accuracy</span>
    <span class="check-detail">Weighted $L^2$ error is below $10^{-2}$ at both checkpoints. The mismatched non-uniform pair outperforms an evenly split uniform grid with the same total number of points.</span>
  </li>
</ul>
~~~

This scenario exercises what multi-domain non-uniform grids are for: different subdomains carry different local feature scales, and the discretization must connect them without introducing artifacts at the seam. The velocity field must remain continuous across the interface; a genuine velocity jump destabilizes the cross-seam WENO stencil continuation. That is a limitation of the current non-conservative, node-centered formulation, not of the grid machinery itself.

<hr class="section-rule">

## Conclusion

The non-uniform WENO foundation in MethodOfLines.jl is now mathematically consistent and operationally complete: dynamic smoothness indicators on stretched stencils, boundary and interface handling, and end-to-end discretization on arbitrary rectilinear grids. The two scenarios above show where that matters. On a uniform grid, under-resolving a steep layer is a failure regime, not a convergence problem. On a clustered grid, the same operator resolves the feature accurately with far fewer global degrees of freedom. Across a multi-domain interface, independently designed non-uniform grids connect seamlessly while tracking a compressing front against an exact solution.

~~~
<div class="cta-box">
  <p class="cta-label">Full implementation</p>
  <p>The mathematical infrastructure, implementation details, and quantitative showcase are collected in <a href="https://github.com/SciML/MethodOfLines.jl/pull/616">PR #616</a>.</p>
</div>
~~~

~~~
<p class="ack">Thank you to my mentor Chris Rackauckas for the technical direction throughout the project.</p>
~~~

~~~
<nav class="series-nav" aria-label="Related posts">
  <p class="series-nav-title">GSoC 2026 series</p>
  <ul>
    <li><a href="/kickoff/">Kickoff</a></li>
    <li><a href="/weno-smoothness/">WENO Smoothness Indicators</a></li>
    <li><a href="/nonuniform-weno-boundaries/">Non-Uniform Boundaries</a></li>
    <li><span class="series-current">Showcase (this post)</span></li>
  </ul>
</nav>
~~~

@@
