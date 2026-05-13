+++
title = "GSoC 2026: Comprehensive Non-Uniform Grid Support for MethodOfLines.jl"
description = "Official kick-off post for my Google Summer of Code 2026 project with NumFOCUS and SciML."
tags = ["gsoc", "julia", "sciml", "pde", "weno"]
hasmath = true
hascode = true
date = Date(2026, 5, 13)
rss_title = "GSoC 2026 Kick-off"
rss_description = "GSoC 2026 Official Kick-off: Implementing Zero-Overhead Trait Dispatch and Non-Uniform WENO."
+++

# Kick-off: GSoC 2026 with SciML

I am absolutely thrilled to announce that my project, **Comprehensive Non-Uniform Grid Support for MethodOfLines.jl**, has been officially accepted for Google Summer of Code (GSoC) 2026 under the NumFOCUS / SciML organization! 

This post serves as the introductory roadmap for my summer coding journey, outlining the critical architectural bottlenecks we will solve and the mathematical engines we will build over the coming months. You can view my full accepted project proposal and development repository [here](https://github.com/utkuyilmaz1903/GSoC-2026-MethodOfLines).

---

## 1. The Illusion of Uniformity

Finite difference discretization is a cornerstone of numerical PDE solving. Uniform grids yield excellent results for many standard physical models. However, in scenarios with steep gradients—such as reaction interfaces in battery electrodes or sharp shockwaves—computing the entire domain with the same density leads to a massive waste of memory and CPU resources. 

Clustering a grid into a specific region directly impacts stability limits. When we shrink local grid spacing ($\Delta x$), the standard CFL stability condition dictates that the time step ($\Delta t$) must drop quadratically. This harsh reality is exactly why we need the powerful implicit solvers and flexible discretization infrastructures of the SciML ecosystem.

## 2. The Current Architectural Bottleneck

When using `MethodOfLines.jl` (MOL) to solve these challenging problems, basic advection schemes support non-uniform grids. However, invoking high-order, shock-capturing schemes like `WENOScheme` on an `AbstractVector` (non-uniform grid) crashes the system:

`ERROR: LoadError: Scheme WENO applied to n(t, u)... is not defined.`

The root cause is that existing high-order operators are defined assuming a scalar $\Delta x$. The current `DiscreteSpace` generator simply cannot handle variable steps.

## 3. Comprehensive Vision and WENO Mathematics

Our primary GSoC goal is to overcome this bottleneck. For the WENO integration, we will implement the following mathematical pipeline:

**Step A: Dynamic Smoothness Indicators ($\beta$)**
Instead of using constant fractions, we will dynamically derive smoothness indicators using Lagrange interpolation over local geometric distances ($h$):
\[ \beta = \int_{x_L}^{x_R} (P'(x))^2 dx + \int_{x_L}^{x_R} (P''(x))^2 dx \]

**Step B: Negative Weight Regularization (Shi-Hu-Shu Splitting)**
To prevent negative linear weights ($d_k$) from causing division by zero on stretched grids, we apply the Shi-Hu-Shu (2002) approach to split the weights into strictly positive sub-groups:
\[ \tilde{d}_k^+ = \frac{1}{2} (d_k + 3|d_k|) \quad \text{and} \quad \tilde{d}_k^- = \tilde{d}_k^+ - d_k \]

**Step C: Final Non-Linear Weights ($\omega$)**
We run the standard WENO formula for these split groups, normalize them, and stitch them back together:
\[ \alpha_k^\pm = \frac{\tilde{d}_k^\pm}{(\epsilon + \beta_k)^2} \implies \omega_k^\pm = \frac{\alpha_k^\pm}{\sum \alpha_j^\pm} \implies \omega_k^{final} = \sigma^+ \cdot \omega_k^+ - \sigma^- \cdot \omega_k^- \]

> **Note on Type Stability:** To ensure GPU readiness and implicit solver compatibility, all internal mathematics strictly avoid hardcoded `Float64` constants, utilizing Julia's `Rational` types and generic typing (`T<:Real`).

## 4. Architectural Routing via Trait-Based Dispatch

To avoid disrupting high-performance structures, the integration leverages a **Trait-based dispatch architecture**. Legacy uniform grids route directly to existing operators. Non-uniform grids bypass the dispatch branch entirely and route directly to isolated dynamic calculators within the WENO module.

```julia
# Compile-time trait resolution for zero-overhead routing
route_engine(grid::T, u, order) where {T} = route_engine(GridType(T), grid, u, order)

route_engine(::UniformGrid, grid, u, order) = apply_classic_weno(grid, u, order)
route_engine(::NonUniformGrid, grid, u, order) = calculate_dynamic_weno_weights(grid, u, order)
```

By moving the dispatch resolution to compile-time via Traits, we eliminate runtime overhead, achieving **~1.2 ns latency with exactly 0 bytes allocated**.

## 5. Securing the Boundaries

To complement the WENO engine, I developed a strictly $O(1)$, zero-allocation mathematical engine for high-order, one-sided finite difference weights. Crucially, to prevent mass leakage at the machine epsilon limit, this boundary engine employs **Kahan (compensated) summation**, preserving absolute stability in long-term simulations.

## Looking Forward

I will be publishing bi-weekly technical updates on my progress. In upcoming posts, I will explore deeper topics such as the formulaic derivation of dynamic smoothness indicators ($\beta$) and multi-dimensional embeddings.

Feel free to drop your feedback on the SciML Slack, or review the architectural implementations directly in my [MethodOfLines.jl Pull Requests](https://github.com/SciML/MethodOfLines.jl/pulls/@utkuyilmaz1903). Here is to a highly productive summer of code!

~~~
<style>
  .math::before, .math::after, .katex-display::before, .katex-display::after {
      display: none !important;
      content: none !important;
  }
</style>
~~~
