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

This post serves as the introductory roadmap for my summer coding journey. You can view my full project proposal [here](https://github.com/utkuyilmaz1903/GSoC-2026-MethodOfLines).

---

## 1. The Illusion of Uniformity

Finite difference discretization is a cornerstone of numerical PDE solving. Uniform grids yield excellent results for many models, but scenarios with steep gradients lead to a massive waste of resources. 

Clustering a grid impacts stability limits. When we shrink local grid spacing ($\Delta x$), the standard CFL condition dictates that the time step ($\Delta t$) must drop quadratically:
\[ \Delta t \le \frac{(\Delta x)^2}{2\alpha} \]

## 2. The Current Architectural Bottleneck

Invoking high-order schemes like `WENOScheme` on an `AbstractVector` (non-uniform grid) currently crashes `MethodOfLines.jl`:
`ERROR: LoadError: Scheme WENO applied to n(t, u)... is not defined.`

## 3. Comprehensive Vision and WENO Mathematics

We will implement a mathematical pipeline to handle dynamic smoothness indicators ($\beta$) and negative weight regularization (Shi-Hu-Shu splitting) to maintain stability on stretched grids.

## 4. Architectural Routing via Trait-Based Dispatch

We leverage a **Trait-based dispatch architecture** to ensure zero-overhead. Initial benchmarks show **~1.2 ns latency with 0 bytes allocated**.

## 5. Securing the Boundaries

To prevent mass leakage, I developed a zero-allocation engine for boundary stencils using **Kahan (compensated) summation**.

## Looking Forward

I will publish bi-weekly updates on my progress. Feel free to drop feedback on the SciML Slack!

~~~
<style>
  .math::before, .math::after, .katex-display::before, .katex-display::after {
      display: none !important;
      content: none !important;
  }
</style>
~~~
