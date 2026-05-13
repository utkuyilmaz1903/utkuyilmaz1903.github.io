+++
title = "GSoC 2026: Comprehensive Non-Uniform Grid Support for MethodOfLines.jl"
description = "Introductory roadmap for my Google Summer of Code 2026 project with SciML."
tags = ["gsoc", "julia", "sciml", "pde", "weno"]
hasmath = true
hascode = true
date = Date(2026, 5, 13)
rss_title = "GSoC 2026 Kick-off"
rss_description = "Implementing Zero-Overhead Trait Dispatch and Non-Uniform WENO."
+++

# Hello everyone!

I am **Utku Yılmaz**, a second-year Computer Engineering student at Muğla 
Sıtkı Koçman University. I am super excited to announce that I have been 
accepted into **Google Summer of Code (GSoC) 2026**! This summer, I will 
be working with the **NumFOCUS / SciML** organization on the 
project: **Comprehensive Non-Uniform Grid Support for MethodOfLines.jl**.

You can view my full accepted project proposal, supplementary materials, 
and development repository [here](https://github.com/utkuyilmaz1903/GSoC-2026-MethodOfLines).

---

## The Engineering Necessity

Finite difference discretization is a cornerstone of numerical PDE solving. 
Building the mathematical infrastructure for non-uniform grids is not just a 
coding task; it is an engineering necessity for simulating real-world physics 
accurately. Uniform grids yield excellent results for many models, but in 
scenarios with steep gradients—such as reaction interfaces in battery electrodes 
or sharp shockwaves—computing the entire domain with the same density leads to 
a massive waste of memory and CPU resources. 

While `MethodOfLines.jl` provides foundational support for non-uniform grids 
via Fornberg weights, high-order schemes like `WENOScheme` currently lack the 
mathematical operators and dispatch logic required for variable spacing. 
Attempting to apply these schemes on clustered domains leads to a hard crash:
`ERROR: LoadError: Scheme WENO applied to n(t, u)... is not defined.`

This summer, I will bridge this architectural gap to enable high-precision 
physical simulations on irregular geometries.

## The Roadmap: Engineering the Discretization Pipeline

To deliver comprehensive, end-to-end non-uniform grid support, my project is 
structured into distinct phases, focusing on both mathematical rigor and 
high-performance Julia code architecture.

### Phase 1: Core Mathematics and Multidimensional Integration (May 25 - July 5)
* **Refactoring Basic Schemes:** We will begin by auditing and refactoring the 
  core advection schemes to natively accept `AbstractVector` grid types, ensuring 
  they handle variable step sizes flawlessly.
* **Isolated WENO Mathematics:** Instead of modifying the core `fornberg.jl`, we 
  will build a dedicated, isolated non-uniform weight calculator entirely within 
  the WENO module. This includes mathematically deriving dynamic smoothness 
  indicators ($\beta$) via Lagrange interpolation and implementing Shi-Hu-Shu 
  (2002) negative weight regularization.
* **Zero-Overhead Dispatch Hub:** I will integrate a Trait-based multiple 
  dispatch infrastructure. This routes legacy uniform grids to existing high-
  performance operators, while automatically redirecting non-uniform arrays to 
  the new isolated engines with absolute zero runtime overhead.
* **Multidimensional Extension:** Finally, the 1D non-uniform WENO stencils 
  will be plugged into the existing multidimensional codegen machinery, backed 
  by rigorous convergence tests to maintain theoretical accuracy.

### Phase 2: Advanced Integration and Capstone (July 6 - August 2)
* **Optimization & Edge Cases:** This phase will focus on optimizing dynamic 
  weight generation and resolving extreme clustering edge cases (e.g., 
  1:1,000,000,000 grid stretching ratios) using Kahan-compensated boundary stencils.
* **Rigorous Benchmarking:** Utilizing `BenchmarkTools.jl`, we will profile 
  execution times against uniform counterparts and formally validate $L_2$ error 
  norms to ensure mathematical correctness.
* **The Capstone Test Case:** To prove the new infrastructure, I will build a 
  multi-domain PDE test case featuring high boundary gradients, inspired by 
  battery interface dynamics. 

### Final Wrap-up (August 3 - August 16)
* The final weeks will be dedicated to cleaning up open PRs, merging the 
  codebase, and expanding the official `MethodOfLines.jl` documentation with 
  brand new tutorials for custom non-uniform grids.

## Upcoming Technical Blog Series

As part of the GSoC journey, I will be publishing bi-weekly technical updates 
on this blog to document my progress. Based on my proposal timeline, here are 
the technical deep-dives you can expect over the summer:

1. **Deriving Non-Uniform Smoothness Indicators ($k$):** Breaking down the math 
   behind replacing classic uniform fractions (like 13/12) with dynamic 
   Lagrange-based indicators.
2. **Dynamic Weight Calculation with Fornberg:** Exploring how to calculate 
   local weights node-by-node without incurring massive memory allocations.
3. **Handling Non-Uniform Boundaries in Multi-Dimensional PDEs:** A look at how 
   univariate spatial derivatives act as drop-in stencil providers for 
   higher-dimensional discretizations.
4. **Profiling Dynamic Stencil Allocations in Julia:** A performance analysis 
   showcasing our zero-overhead Trait dispatch system.
5. **$L_2$ Error Norm Convergence Analysis:** Formal mathematical validation 
   proving that our non-uniform WENO scheme maintains its theoretical $O(\Delta x^3)$ 
   convergence order.

I would like to give a huge thanks to my mentors, **Chris Rackauckas** and 
**Alex Jones**, for their guidance and this amazing opportunity to contribute 
to the bleeding edge of scientific machine learning. 

Feel free to drop your feedback on the SciML Slack, or track my progress directly 
on GitHub. Let's build something great!

~~~
<style>
  .math::before, .math::after, .katex-display::before, .katex-display::after {
      display: none !important;
      content: none !important;
  }
</style>
~~~
